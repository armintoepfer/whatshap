#include <stdexcept>
#include <cassert>
#include <limits>
#include <fstream>
#include <array>
#include <algorithm>
#include <cmath>
#include <vector>
#include <cstdint>

#include "genotype.h"
#include "pedigreecolumncostcomputer.h"
#include "pedigreedptable.h"

using namespace std;

namespace {

/** The columns of a DP for one individual (see PedigreeDPTable::compute_table_single_individual).
 *  Entry j of a column is its j-th active read (read id order), bit j of a bipartition index. */
struct single_individual_columns_t {
	// per column: the number of active reads, the reads shared with the previous column,
	// and the offset of its entries
	vector<unsigned int> size;
	vector<unsigned int> backward_width;
	vector<size_t> offset;
	// per entry: the phred score of a REF allele (0 otherwise), of an ALT allele (0 otherwise),
	// and its bit in the forward projection (-1: the read ends in this column)
	vector<unsigned int> ref_score;
	vector<unsigned int> alt_score;
	vector<int> forward_bit;
	unsigned int max_size = 0;
};

/** The step of the Gray code of GrayCodes that visits a bipartition index (the inverse Gray code). */
inline unsigned int gray_rank(unsigned int index) {
	index ^= index >> 1;
	index ^= index >> 2;
	index ^= index >> 4;
	index ^= index >> 8;
	index ^= index >> 16;
	return index;
}

/** The index of the lowest set bit of x > 0: the bit that GrayCodes changes in step x. */
inline int lowest_bit(unsigned int x) {
#if defined(__GNUC__)
	return __builtin_ctz(x);
#else
	int i = 0;
	while ((x & 1u) == 0) {
		x >>= 1;
		++i;
	}
	return i;
#endif
}

/** One column of the DP in the order of ColumnIndexingIterator (Gray code), for columns of any size.
 *  previous: the projection of the previous column; projection, backtrace: the forward projection
 *  of this column and the index that gives each entry (not written for the last column).
 *  For the last column, returns the first index with the minimum value in optimal_index. */
template <class Index>
void single_individual_column_gray(const single_individual_columns_t& columns, size_t c, bool last, const vector<unsigned int>& previous, vector<unsigned int>& projection, Index* backtrace, unsigned int& optimal_score, unsigned int& optimal_index) {
	const unsigned int n = columns.size[c];
	const size_t o = columns.offset[c];
	const unsigned int backward_mask = (((unsigned int)1) << columns.backward_width[c]) - 1;
	// the costs of the entries in partition 0 and 1, by allele
	unsigned int ref0 = 0, alt0 = 0, ref1 = 0, alt1 = 0;
	for (unsigned int j = 0; j < n; ++j) {
		ref0 += columns.ref_score[o + j];
		alt0 += columns.alt_score[o + j];
	}
	unsigned int forward = 0;
	const unsigned int states = ((unsigned int)1) << n;
	for (unsigned int k = 0; k < states; ++k) {
		const unsigned int index = k ^ (k >> 1);
		if (k > 0) {
			const int b = lowest_bit(k);
			if (((index >> b) & 1) != 0) {
				ref0 -= columns.ref_score[o + b];
				alt0 -= columns.alt_score[o + b];
				ref1 += columns.ref_score[o + b];
				alt1 += columns.alt_score[o + b];
			} else {
				ref1 -= columns.ref_score[o + b];
				alt1 -= columns.alt_score[o + b];
				ref0 += columns.ref_score[o + b];
				alt0 += columns.alt_score[o + b];
			}
			if (columns.forward_bit[o + b] >= 0) {
				forward ^= ((unsigned int)1) << columns.forward_bit[o + b];
			}
		}
		// the two allele assignments: partition 0 ALT and partition 1 REF, or the opposite
		unsigned int value = std::min(ref0 + alt1, alt0 + ref1);
		if (c > 0) {
			const unsigned int previous_value = previous[index & backward_mask];
			value = (previous_value == numeric_limits<unsigned int>::max()) ? previous_value : value + previous_value;
		}
		if (last) {
			if (value < optimal_score) {
				optimal_score = value;
				optimal_index = index;
			}
		} else if (value < projection[forward]) {
			projection[forward] = value;
			backtrace[forward] = (Index)index;
		}
	}
}

/** The same column as array operations, for columns of at most 16 reads. The reads shared with the
 *  previous column are the low bits of an index, so the value of each index is its cost plus the
 *  previous projection at its low bits. The forward projection then takes, for each read that ends,
 *  the better index of each pair that differs in its bit. The Gray code order keeps the first
 *  minimum (strict <), so a tie goes to the index of the smaller Gray rank, here too. */
template <class Index>
void single_individual_column_array(const single_individual_columns_t& columns, size_t c, bool last, const vector<unsigned int>& previous, vector<unsigned int>& projection, Index* backtrace, unsigned int& optimal_score, unsigned int& optimal_index, vector<int>& diff, vector<unsigned int>& value, vector<unsigned int>& origin) {
	const unsigned int n = columns.size[c];
	const size_t o = columns.offset[c];
	const size_t states = ((size_t)1) << n;
	const size_t shared = ((size_t)1) << columns.backward_width[c];
	// cost = min(ref0 + alt1, alt0 + ref1) = min(ref_total - d, alt_total + d), d = ref1 - alt1
	long long ref_total = 0, alt_total = 0;
	diff.assign(states, 0);
	for (unsigned int j = 0; j < n; ++j) {
		ref_total += columns.ref_score[o + j];
		alt_total += columns.alt_score[o + j];
		const int d = (int)columns.ref_score[o + j] - (int)columns.alt_score[o + j];
		const size_t half = ((size_t)1) << j;
		for (size_t i = 0; i < half; ++i) {
			diff[half + i] = diff[i] + d;
		}
	}
	value.resize(states);
	origin.resize(states);
	for (size_t high = 0; high < states; high += shared) {
		for (size_t low = 0; low < shared; ++low) {
			const long long d = diff[high + low];
			const unsigned int cost = (unsigned int)std::min(ref_total - d, alt_total + d);
			unsigned int v = cost;
			if (c > 0) {
				v = (previous[low] == numeric_limits<unsigned int>::max()) ? previous[low] : cost + previous[low];
			}
			value[high + low] = v;
			origin[high + low] = (unsigned int)(high + low);
		}
	}
	if (last) {
		size_t best = 0;
		for (size_t i = 1; i < states; ++i) {
			if ((value[i] < value[best]) || ((value[i] == value[best]) && (gray_rank(origin[i]) < gray_rank(origin[best])))) {
				best = i;
			}
		}
		if (value[best] < optimal_score) {
			optimal_score = value[best];
			optimal_index = origin[best];
		}
		return;
	}
	// remove the bit of each read that ends here, the highest first (the lower bits keep their places)
	size_t size = states;
	for (unsigned int j = n; j-- > 0;) {
		if (columns.forward_bit[o + j] >= 0) {
			continue;
		}
		const size_t bit = ((size_t)1) << j;
		size_t w = 0;
		for (size_t block = 0; block < size; block += 2 * bit) {
			for (size_t i = block; i < block + bit; ++i) {
				const bool second = (value[i + bit] < value[i]) || ((value[i + bit] == value[i]) && (gray_rank(origin[i + bit]) < gray_rank(origin[i])));
				value[w] = second ? value[i + bit] : value[i];
				origin[w] = second ? origin[i + bit] : origin[i];
				++w;
			}
		}
		size = w;
	}
	assert(size == projection.size());
	for (size_t i = 0; i < size; ++i) {
		projection[i] = value[i];
		backtrace[i] = (Index)origin[i];
	}
}

/** Computes the columns [from, to) of the DP for one individual. previous holds the projection of
 *  column from-1 (empty for column 0) and then the projection of column to-1. backtrace and offset
 *  get the backtrace of each column c < to (except the last column of the table), at offset[c - from].
 *  The last column of the table sets optimal_score and optimal_index. */
template <class Index>
void single_individual_forward(const single_individual_columns_t& columns, size_t from, size_t to, vector<unsigned int>& previous, vector<Index>& backtrace, vector<size_t>& offset, unsigned int& optimal_score, unsigned int& optimal_index, vector<unsigned int>& projection, vector<int>& diff, vector<unsigned int>& value, vector<unsigned int>& origin) {
	const size_t column_count = columns.size.size();
	offset.assign(to - from, 0);
	size_t total = 0;
	for (size_t c = from; (c < to) && (c + 1 < column_count); ++c) {
		offset[c - from] = total;
		total += ((size_t)1) << columns.backward_width[c + 1];
	}
	backtrace.resize(total);
	for (size_t c = from; c < to; ++c) {
		const bool last = (c + 1 == column_count);
		Index* column_backtrace = nullptr;
		if (!last) {
			projection.assign(((size_t)1) << columns.backward_width[c + 1], numeric_limits<unsigned int>::max());
			column_backtrace = backtrace.data() + offset[c - from];
		}
		if (columns.size[c] <= 16) {
			single_individual_column_array<Index>(columns, c, last, previous, projection, column_backtrace, optimal_score, optimal_index, diff, value, origin);
		} else {
			single_individual_column_gray<Index>(columns, c, last, previous, projection, column_backtrace, optimal_score, optimal_index);
		}
		if (!last) {
			previous.swap(projection);
		}
	}
}

/** The forward pass and the backtrace of the DP for one individual. Fills index_path with the
 *  bipartition index of each column and returns the optimal score. As compute_table, it keeps the
 *  projection of every k-th column (k = sqrt(#columns)) and computes the backtrace of a block of
 *  k columns again when the backtrace needs it. */
template <class Index>
unsigned int solve_single_individual(const single_individual_columns_t& columns, vector<index_and_inheritance_t>& index_path) {
	const size_t column_count = columns.size.size();
	const size_t k = std::max((size_t)1, (size_t)sqrt(column_count));
	const size_t block_count = (column_count + k - 1) / k;
	// checkpoints[b]: the projection that enters column b*k
	vector<vector<unsigned int> > checkpoints(block_count);
	vector<unsigned int> previous;
	vector<Index> backtrace;
	vector<size_t> offset;
	vector<unsigned int> projection, value, origin;
	vector<int> diff;
	unsigned int optimal_score = numeric_limits<unsigned int>::max();
	unsigned int optimal_index = 0;
	for (size_t b = 0; b < block_count; ++b) {
		checkpoints[b] = previous;
		single_individual_forward<Index>(columns, b * k, std::min((b + 1) * k, column_count), previous, backtrace, offset, optimal_score, optimal_index, projection, diff, value, origin);
	}

	// backtrace and offset hold the last block now
	index_path.assign(column_count, index_and_inheritance_t());
	unsigned int index = optimal_index;
	size_t block = block_count - 1;
	for (size_t c = column_count; c-- > 0;) {
		index_path[c].index = index;
		if (c == 0) {
			break;
		}
		if (c - 1 < block * k) {
			// column c-1 is in the previous block: compute its backtrace again
			--block;
			previous = checkpoints[block];
			unsigned int unused_score = numeric_limits<unsigned int>::max();
			unsigned int unused_index = 0;
			single_individual_forward<Index>(columns, block * k, (block + 1) * k, previous, backtrace, offset, unused_score, unused_index, projection, diff, value, origin);
		}
		const unsigned int backward_mask = (((unsigned int)1) << columns.backward_width[c]) - 1;
		index = backtrace[offset[c - 1 - block * k] + (index & backward_mask)];
	}
	return optimal_score;
}

/** One column of PedigreeDPTable::compute_column_fast: its sizes and its tables (rows of
 *  "configurations" values). previous is nullptr for the first column; projection and the two
 *  backtraces are nullptr for the last column. */
struct pedigree_column_t {
	size_t configurations;
	size_t entries;
	size_t mask_count;
	size_t width;
	unsigned int backward_width;
	const unsigned int* previous;
	unsigned int* projection;
	unsigned int* index_backtrace;
	unsigned int* transmission_backtrace;
};

/** The optimum of the last column, as in PedigreeDPTable. */
struct pedigree_optimum_t {
	unsigned int score;
	unsigned int index;
	unsigned int transmission_value;
	unsigned int previous_transmission_value;
};

/** The bipartitions of one column for PedigreeDPTable::compute_column_fast, in the order of
 *  ColumnIndexingIterator (Gray code). T is the number of transmission configurations, or 0 when
 *  it is only known at run time. value: the cost of each mask for bipartition 0 (and a last value of
 *  infinity); delta: the change of each mask cost when entry j moves to haplotype 0 or 1;
 *  assignment: the allowed masks of each configuration (column.width each); recombination: T x T. */
template <size_t T>
void pedigree_column(const pedigree_column_t& column, vector<unsigned int>& value, const vector<unsigned int>& delta, const vector<unsigned int>& assignment, const vector<unsigned int>& recombination, const vector<unsigned int>& forward_flip, pedigree_optimum_t& optimum) {
	const size_t configurations = (T > 0) ? T : column.configurations;
	const size_t mask_count = column.mask_count;
	const size_t width = column.width;
	const unsigned int infinity = numeric_limits<unsigned int>::max();

	// For each entry b of the previous projection and each configuration i: the minimum over the
	// configurations j of the previous cost plus the recombination cost from j to i, and the first j
	// that gives it. The cost of i in this column is the same for all j, so compute_column_reference
	// finds the same minimum and the same j for each bipartition (all costs are far below 2^32).
	const size_t previous_size = (column.previous != nullptr) ? (((size_t)1) << column.backward_width) : 1;
	vector<unsigned int> previous_min(previous_size * configurations);
	vector<unsigned int> previous_min_index(previous_size * configurations);
	const vector<unsigned int> no_previous_cost(configurations, 0);
	for (size_t b = 0; b < previous_size; ++b) {
		const unsigned int* previous_cost = (column.previous != nullptr) ? column.previous + b * configurations : no_previous_cost.data();
		for (size_t i = 0; i < configurations; ++i) {
			const unsigned int* recombination_to_i = &recombination[i * configurations];
			unsigned int min = infinity;
			unsigned int min_index = 0;
			for (size_t j = 0; j < configurations; ++j) {
				const unsigned int val = (previous_cost[j] < infinity) ? previous_cost[j] + recombination_to_i[j] : infinity;
				const bool better = val < min;
				min = better ? val : min;
				min_index = better ? (unsigned int)j : min_index;
			}
			previous_min[b * configurations + i] = min;
			previous_min_index[b * configurations + i] = min_index;
		}
	}

	const unsigned int backward_mask = (column.previous != nullptr) ? (((unsigned int)1) << column.backward_width) - 1 : 0;
	vector<unsigned int> min_value(configurations);
	unsigned int forward_index = 0;
	const unsigned int states = ((unsigned int)1) << column.entries;
	for (unsigned int k = 0; k < states; ++k) {
		const unsigned int index = k ^ (k >> 1);
		if (k > 0) {
			// entry b moves to haplotype (index >> b) & 1 of its individual
			const int b = lowest_bit(k);
			const unsigned int* d = &delta[(2 * b + ((index >> b) & 1)) * mask_count];
			for (size_t m = 0; m < mask_count; ++m) {
				value[m] += d[m];
			}
			forward_index ^= forward_flip[b];
		}

		const size_t previous_offset = (index & backward_mask) * configurations;
		const unsigned int* previous = &previous_min[previous_offset];
		unsigned int min_current_cost = infinity;
		for (size_t i = 0; i < configurations; ++i) {
			unsigned int current_cost = infinity;
			for (size_t a = 0; a < width; ++a) {
				current_cost = std::min(current_cost, value[assignment[i * width + a]]);
			}
			min_current_cost = std::min(min_current_cost, current_cost);
			min_value[i] = ((current_cost < infinity) && (previous[i] < infinity)) ? current_cost + previous[i] : infinity;
		}
		if (min_current_cost == infinity) {
			throw std::runtime_error("Error: Mendelian conflict");
		}

		const unsigned int* min_recomb_index = &previous_min_index[previous_offset];
		if (column.projection == nullptr) {
			for (size_t i = 0; i < configurations; ++i) {
				if (min_value[i] < optimum.score) {
					optimum.score = min_value[i];
					optimum.index = index;
					optimum.transmission_value = i;
					optimum.previous_transmission_value = min_recomb_index[i];
				}
			}
		} else {
			unsigned int* projection = column.projection + forward_index * configurations;
			unsigned int* index_backtrace = column.index_backtrace + forward_index * configurations;
			unsigned int* transmission_backtrace = column.transmission_backtrace + forward_index * configurations;
			for (size_t i = 0; i < configurations; ++i) {
				const bool better = min_value[i] < projection[i];
				projection[i] = better ? min_value[i] : projection[i];
				index_backtrace[i] = better ? index : index_backtrace[i];
				transmission_backtrace[i] = better ? min_recomb_index[i] : transmission_backtrace[i];
			}
		}
	}
}

} // namespace

PedigreeDPTable::PedigreeDPTable(ReadSet* read_set, const vector<unsigned int>& recombcost, const Pedigree* pedigree, bool distrust_genotypes, const vector<unsigned int>* positions, bool allow_fast_dp) :
	read_set(read_set),
	recombcost(recombcost),
	pedigree(pedigree),
	distrust_genotypes(distrust_genotypes),
	allow_fast_dp(allow_fast_dp),
	optimal_score(0u),
	optimal_score_index(0u),
	input_column_iterator(*read_set, positions)
{
	read_set->reassignReadIds();

	// create all pedigree partitions
	for (size_t i=0; i<std::pow(4, pedigree->triple_count()); ++i) {
		pedigree_partitions.push_back(new PedigreePartitions(*pedigree, i));
	}

	// translate all individual ids to individual indices
	for (size_t i=0; i<read_set->size(); ++i) {
		read_sources.push_back(pedigree->id_to_index(read_set->get(i)->getSampleID()));
	}

	compute_table();
}


PedigreeDPTable::~PedigreeDPTable() {
	init(projection_column_table, 0);
	init(index_backtrace_table, 0);
	init(transmission_backtrace_table, 0);
	init(indexers, 0);
	init(pedigree_partitions, 0);
}


unique_ptr<vector<unsigned int> > PedigreeDPTable::extract_read_ids(const vector<const Entry *>& entries) {
	unique_ptr<vector<unsigned int> > read_ids(new vector<unsigned int>());
	for (size_t i=0; i<entries.size(); ++i) {
		read_ids->push_back(entries[i]->get_read_id());
	}
	return read_ids;
}


size_t PedigreeDPTable::popcount(size_t x) {
	unsigned int count = 0;
	for (;x; x >>= 1) {
		count += x & 1;
	}
	return count;
}


void PedigreeDPTable::clear_table() {
	size_t column_count = input_column_iterator.get_column_count();

	init(projection_column_table, column_count);
	init(index_backtrace_table, column_count);
	init(transmission_backtrace_table, column_count);
	init(indexers, column_count);

	index_path.clear();

	optimal_score = numeric_limits<unsigned int>::max();
	optimal_score_index = 0;
	optimal_transmission_value = 0;
	previous_transmission_value = 0;
}


void PedigreeDPTable::compute_table() {
	clear_table();

	// empty read-set, nothing to phase, so MEC score is 0
	if (input_column_iterator.get_column_count() == 0) {
		optimal_score = 0;
		optimal_score_index = 0;
		return;
	}

	if (allow_fast_dp && is_single_individual_het()) {
		compute_table_single_individual();
		return;
	}

	input_column_iterator.jump_to_column(0);
	unique_ptr<vector<const Entry *> > current_input_column;
	unique_ptr<vector<const Entry *> > next_input_column;
	// get the next column ahead of time
	next_input_column = input_column_iterator.get_next();
	unique_ptr<vector<unsigned int> > next_read_ids = extract_read_ids(*next_input_column);
	ColumnIndexingScheme* next_indexer = new ColumnIndexingScheme(0, *next_read_ids);
	indexers[0] = next_indexer;

	// forward pass: create a sparse table, storing values at every sqrt(#columns)-th position,
	size_t k = (size_t)sqrt(input_column_iterator.get_column_count());
	for (size_t column_index=0; column_index<input_column_iterator.get_column_count(); ++column_index) {
		// make former next column the current one
		current_input_column = std::move(next_input_column);
		unique_ptr<vector<unsigned int> > current_read_ids = std::move(next_read_ids);

		ColumnIndexingScheme* current_indexer = next_indexer;
		// peek ahead and get the next column
		if (input_column_iterator.has_next()) {
			next_input_column = input_column_iterator.get_next();
			next_read_ids = extract_read_ids(*next_input_column);
			next_indexer = new ColumnIndexingScheme(current_indexer,*next_read_ids);
			current_indexer->set_next_column(next_indexer);
			indexers[column_index + 1] = next_indexer;
		} else {
			assert(next_input_column.get() == 0);
			assert(next_read_ids.get() == 0);
			next_indexer = 0;
		}

		compute_column(column_index, std::move(current_input_column));

		// determine whether to delete previous column (to save space)
		if ((k>1) && (column_index > 0) && (((column_index-1)%k) != 0)) {
			delete index_backtrace_table[column_index-1];
			delete transmission_backtrace_table[column_index-1];
			delete projection_column_table[column_index-1];
			index_backtrace_table[column_index-1] = nullptr;
			transmission_backtrace_table[column_index-1] = nullptr;
			projection_column_table[column_index-1] = nullptr;
		}
	}

	// perform a backtrace to get optimal path
	index_path.assign(indexers.size(), index_and_inheritance_t());
	index_and_inheritance_t v;
	unsigned int prev_inheritance_value = previous_transmission_value;
	v.index = optimal_score_index;
	v.inheritance_value = optimal_transmission_value;
	index_path[indexers.size()-1] = v;
	for(size_t i = indexers.size()-1; i > 0; --i) { // backtrack through table
		// ensure that index_backtrace_table[i-1] and transmission_backtrace_table[i-1] exist
		if (projection_column_table[i-1] == nullptr) {
			// compute index of last previous column that has been stored
			size_t j = (i-1) / k * k;
			assert(projection_column_table[j] != nullptr);
			for (j=j+1; j<i; ++j) {
				compute_column(j);
			}
		}
		// compute index and transmission value for the current column
		unique_ptr<ColumnIndexingIterator> iterator = indexers[i]->get_iterator();
		unsigned int backtrace_index = iterator->index_backward_projection(v.index);
		v.index = index_backtrace_table[i-1]->at(backtrace_index, prev_inheritance_value);
		v.inheritance_value = prev_inheritance_value;
		prev_inheritance_value = transmission_backtrace_table[i-1]->at(backtrace_index, v.inheritance_value);
		index_path[i-1] = v;
		// free parts of the DP table no longer needed
		if (i%k == 0) {
			for (size_t j=i; (j<i+k) && (j<input_column_iterator.get_column_count()-1); ++j) {
				assert(projection_column_table[j] != nullptr);
				delete index_backtrace_table[j];
				delete transmission_backtrace_table[j];
				delete projection_column_table[j];
				index_backtrace_table[j] = nullptr;
				transmission_backtrace_table[j] = nullptr;
				projection_column_table[j] = nullptr;
			}
		}
	}
}


bool PedigreeDPTable::is_single_individual_het() {
	if (distrust_genotypes || (pedigree->size() != 1) || (pedigree->triple_count() != 0)) {
		return false;
	}
	const Genotype het(vector<uint32_t>{0, 1});
	for (size_t i = 0; i < input_column_iterator.get_positions()->size(); ++i) {
		const Genotype* genotype = pedigree->get_genotype(0, i);
		if ((genotype == nullptr) || (*genotype != het)) {
			return false;
		}
	}
	return true;
}


void PedigreeDPTable::compute_table_single_individual() {
	const size_t column_count = input_column_iterator.get_column_count();
	single_individual_columns_t columns;
	columns.size.reserve(column_count);
	columns.backward_width.reserve(column_count);
	columns.offset.reserve(column_count);

	// the same indexers as compute_table (get_optimal_partitioning reads them)
	input_column_iterator.jump_to_column(0);
	for (size_t column_index = 0; column_index < column_count; ++column_index) {
		unique_ptr<vector<const Entry *> > column = input_column_iterator.get_next();
		unique_ptr<vector<unsigned int> > read_ids = extract_read_ids(*column);
		ColumnIndexingScheme* previous_indexer = (column_index > 0) ? indexers[column_index - 1] : nullptr;
		indexers[column_index] = new ColumnIndexingScheme(previous_indexer, *read_ids);
		if (previous_indexer != nullptr) {
			previous_indexer->set_next_column(indexers[column_index]);
		}
		columns.size.push_back(column->size());
		columns.backward_width.push_back(indexers[column_index]->get_backward_projection_width());
		columns.offset.push_back(columns.ref_score.size());
		columns.max_size = std::max(columns.max_size, (unsigned int)column->size());
		for (const Entry* entry : *column) {
			switch (entry->get_allele_type()) {
			case Entry::REF_ALLELE:
				columns.ref_score.push_back(entry->get_phred_score());
				columns.alt_score.push_back(0);
				break;
			case Entry::ALT_ALLELE:
				columns.ref_score.push_back(0);
				columns.alt_score.push_back(entry->get_phred_score());
				break;
			case Entry::BLANK:
				columns.ref_score.push_back(0);
				columns.alt_score.push_back(0);
				break;
			default:
				assert(false);
			}
			columns.forward_bit.push_back(-1);
		}
	}
	// the bit of each entry in the forward projection of its column
	for (size_t column_index = 0; column_index + 1 < column_count; ++column_index) {
		const vector<unsigned int>* mask = indexers[column_index]->get_forward_projection_mask();
		for (size_t j = 0; j < mask->size(); ++j) {
			columns.forward_bit[columns.offset[column_index] + j] = (int)mask->at(j);
		}
	}

	// The backtrace takes 16 bits per entry when no column has more than 16 active reads
	// (the read selection keeps a coverage of 15 by default).
	if (columns.max_size <= 16) {
		optimal_score = solve_single_individual<uint16_t>(columns, index_path);
	} else {
		optimal_score = solve_single_individual<uint32_t>(columns, index_path);
	}
	optimal_score_index = index_path.back().index;
	optimal_transmission_value = 0;
	previous_transmission_value = 0;
}


void PedigreeDPTable::compute_column(size_t column_index, unique_ptr<vector<const Entry*>> current_input_column) {
	assert(column_index < input_column_iterator.get_column_count());

	// check whether requested column is already there
	if (projection_column_table[column_index] != nullptr) {
		assert(index_backtrace_table[column_index] != nullptr);
		assert(transmission_backtrace_table[column_index] != nullptr);
		return;
	}

	ColumnIndexingScheme* current_indexer = indexers[column_index];
	assert(current_indexer != nullptr);

	// compute the number of different transmission vectors
	unsigned int transmission_configurations = std::pow(4, pedigree->triple_count());

	// if current input column was not provided, then create it
	if (current_input_column.get() == nullptr) {
		input_column_iterator.jump_to_column(column_index);
		current_input_column = input_column_iterator.get_next();
	}

	// obtain previous projection column (which is assumed to have been already computed)
	Vector2D<unsigned int>* previous_projection_column = nullptr;
	if (column_index > 0) {
		previous_projection_column = projection_column_table[column_index - 1];
	}

	// initialize forward projection column and associated backtrace columns,
	// if existing (i.e. if not last column)
	Vector2D<unsigned int>* current_projection_column = nullptr;
	Vector2D<unsigned int>* transmission_backtrace_column = nullptr;
	Vector2D<unsigned int>* index_backtrace_column = nullptr;
	if (column_index + 1 < input_column_iterator.get_column_count()) {
		current_projection_column = new Vector2D<unsigned int>(
			current_indexer->forward_projection_size(),
			transmission_configurations,
			numeric_limits<unsigned int>::max()
		);
		transmission_backtrace_column = new Vector2D<unsigned int>(
			current_indexer->forward_projection_size(),
			transmission_configurations,
			numeric_limits<unsigned int>::max()
		);
		index_backtrace_column = new Vector2D<unsigned int>(
			current_indexer->forward_projection_size(),
			transmission_configurations,
			numeric_limits<unsigned int>::max()
		);
	}

	if (allow_fast_dp) {
		compute_column_fast(column_index, *current_input_column, previous_projection_column, current_projection_column, index_backtrace_column, transmission_backtrace_column);
	} else {
		compute_column_reference(column_index, *current_input_column, previous_projection_column, current_projection_column, index_backtrace_column, transmission_backtrace_column);
	}

	// if not last column, then store computed tables
	if (current_projection_column != 0) {
		index_backtrace_table[column_index] = index_backtrace_column;
		transmission_backtrace_table[column_index] = transmission_backtrace_column;
		projection_column_table[column_index] = current_projection_column;
	}
}


void PedigreeDPTable::compute_column_reference(size_t column_index, const vector<const Entry*>& current_input_column, const Vector2D<unsigned int>* previous_projection_column, Vector2D<unsigned int>* current_projection_column, Vector2D<unsigned int>* index_backtrace_column, Vector2D<unsigned int>* transmission_backtrace_column) {
	ColumnIndexingScheme* current_indexer = indexers[column_index];
	unsigned int transmission_configurations = pedigree_partitions.size();

	// reserve memory for the current DP column
	Vector2D<unsigned int> dp_column(current_indexer->column_size(), transmission_configurations, 0);

	// create column cost computers
	vector<PedigreeColumnCostComputer> cost_computers;
	cost_computers.reserve(transmission_configurations);
	for(unsigned int i = 0; i < transmission_configurations; ++i) {
		cost_computers.emplace_back(current_input_column, column_index, read_sources, pedigree, *pedigree_partitions[i], distrust_genotypes);
	}

	// iterate over all bipartitions
	unique_ptr<ColumnIndexingIterator> iterator = current_indexer->get_iterator();
	while (iterator->has_next()) {
		int bit_changed = -1;
		iterator->advance(&bit_changed);
		if (bit_changed >= 0) {
			for(auto& cost_computer : cost_computers) {
				cost_computer.update_partitioning(bit_changed);
			}
		} else {
			for(auto& cost_computer : cost_computers) {
				cost_computer.set_partitioning(iterator->get_partition());
			}
		}

		// Determine index in backward projection column from where to fetch the previous cost
		size_t backward_projection_index = 0;
		if (column_index > 0) {
			backward_projection_index = iterator->get_backward_projection();
		}
		// Determine index in the current DP column to be written
		size_t current_index = iterator->get_index();

		// Compute aggregate cost based on cost in previous and cost in current column
		vector<unsigned int> min_recomb_index(transmission_configurations);
		bool found_valid_transmission_vector = false;
		for (size_t i = 0; i < transmission_configurations; ++i) {
			// Compute cost incurred by current cell of DP table
			unsigned int current_cost = cost_computers[i].get_cost();
			unsigned int min = numeric_limits<unsigned int>::max();
			size_t min_index = 0;
			if (current_cost < numeric_limits<unsigned int>::max()) {
				found_valid_transmission_vector = true;
			}
			for (size_t j = 0; j < transmission_configurations; ++j) {
				// Step 1: add up cost from current_cost column and previous columns
				unsigned int val;
				unsigned int previous_cost = 0;
				if (column_index > 0) {
					previous_cost = previous_projection_column->at(backward_projection_index,j);
				}
				if ((current_cost < numeric_limits<unsigned int>::max()) && (previous_cost < numeric_limits<unsigned int>::max())) {
					val = current_cost + previous_cost;
				} else {
					val = numeric_limits<unsigned int>::max();
				}
				// Step 2: add further cost incurred by recombination
				// change in bit 0 --> recombination in mother
				size_t x = i ^ j; // count the number of bits set in x

				if (val < numeric_limits<unsigned int>::max()) {
					val += popcount(x) * recombcost[column_index];
				}

				// check for new minimum
				if (val < min) {
					min = val;
					min_index = j;
				}
			}
			dp_column.set(current_index, i, min);
			min_recomb_index[i] = min_index;
		}
		if (!found_valid_transmission_vector) {
			throw std::runtime_error("Error: Mendelian conflict");
		}

		// if last DP column, then check for new optimal score, otherwise update forward projection and backtrace columns
		if (current_projection_column == 0) {
			// update running optimal score index
			for (size_t i = 0; i < transmission_configurations; ++i) {
				if (dp_column.at(current_index, i) < optimal_score) {
					optimal_score = dp_column.at(current_index, i);
					optimal_score_index = iterator->get_index();
					optimal_transmission_value = i;
					previous_transmission_value = min_recomb_index[i];
				}
			}
		} else {
			unsigned int forward_index = iterator->get_forward_projection();
			unsigned int it_idx = iterator->get_index();
			for (unsigned int i = 0; i < transmission_configurations; ++i) {
				if (dp_column.at(current_index, i) < current_projection_column->at(forward_index,i)) {
					current_projection_column->set(forward_index, i, dp_column.at(current_index, i));
					index_backtrace_column->set(forward_index, i, it_idx);
					transmission_backtrace_column->set(forward_index,i, min_recomb_index[i]);
				}
			}
		}
	}
}


void PedigreeDPTable::compute_column_fast(size_t column_index, const vector<const Entry*>& current_input_column, const Vector2D<unsigned int>* previous_projection_column, Vector2D<unsigned int>* current_projection_column, Vector2D<unsigned int>* index_backtrace_column, Vector2D<unsigned int>* transmission_backtrace_column) {
	ColumnIndexingScheme* current_indexer = indexers[column_index];
	const size_t transmission_configurations = pedigree_partitions.size();
	const size_t n = current_input_column.size();
	assert(n < 32);
	const size_t individuals = pedigree->size();

	// The allele assignments of PedigreeColumnCostComputer, as a mask over the haplotypes of all
	// individuals (bit 2*individual+haplotype). Each mask occurs once, with its genotype cost;
	// assignments[i] lists those that are allowed for transmission configuration i.
	vector<unsigned int> masks;
	vector<unsigned int> genotype_costs;
	vector<vector<unsigned int> > assignments(transmission_configurations);
	for (size_t i = 0; i < transmission_configurations; ++i) {
		const PedigreePartitions& partitions = *pedigree_partitions[i];
		for (unsigned int a = 0; a < (((unsigned int)1) << partitions.count()); ++a) {
			bool genotypes_compatible = true;
			unsigned int cost = 0;
			unsigned int mask = 0;
			for (size_t individuals_index = 0; individuals_index < individuals; ++individuals_index) {
				unsigned int allele0 = (a >> partitions.haplotype_to_partition(individuals_index, 0)) & 1;
				unsigned int allele1 = (a >> partitions.haplotype_to_partition(individuals_index, 1)) & 1;
				mask |= (allele0 << (2 * individuals_index)) | (allele1 << (2 * individuals_index + 1));
				Genotype genotype(vector<unsigned int>{allele0, allele1});
				if (distrust_genotypes) {
					const PhredGenotypeLikelihoods* gls = pedigree->get_genotype_likelihoods(individuals_index, column_index);
					assert(gls != nullptr);
					cost += gls->get(genotype);
				} else {
					const Genotype* true_genotype = pedigree->get_genotype(individuals_index, column_index);
					if (genotype != (*true_genotype)) {
						genotypes_compatible = false;
						break;
					}
				}
			}
			if (!genotypes_compatible) {
				continue;
			}
			const size_t m = std::find(masks.begin(), masks.end(), mask) - masks.begin();
			if (m == masks.size()) {
				masks.push_back(mask);
				genotype_costs.push_back(cost);
			}
			assignments[i].push_back(m);
		}
	}

	// per entry: its individual, and its cost on a haplotype with allele 0 (an ALT entry) or 1 (a REF entry)
	vector<unsigned int> entry_individual(n);
	vector<array<unsigned int, 2> > entry_cost(n, {0, 0});
	for (size_t j = 0; j < n; ++j) {
		const Entry& entry = *current_input_column[j];
		entry_individual[j] = read_sources[entry.get_read_id()];
		switch (entry.get_allele_type()) {
		case Entry::REF_ALLELE:
			entry_cost[j][1] = entry.get_phred_score();
			break;
		case Entry::ALT_ALLELE:
			entry_cost[j][0] = entry.get_phred_score();
			break;
		case Entry::BLANK:
			break;
		default:
			assert(false);
		}
	}

	const unsigned int infinity = numeric_limits<unsigned int>::max();

	// value[m]: the cost of mask m for the current bipartition (bit j: the haplotype of entry j).
	// Bipartition 0 puts all entries on haplotype 0. The sums are modulo 2^32 as in
	// PedigreeColumnCostComputer, so the order of the additions does not change them.
	// The last value (index mask_count) stays infinity: it stands for "no allowed assignment".
	const size_t mask_count = masks.size();
	vector<unsigned int> value(genotype_costs);
	for (size_t m = 0; m < mask_count; ++m) {
		for (size_t j = 0; j < n; ++j) {
			value[m] += entry_cost[j][(masks[m] >> (2 * entry_individual[j])) & 1];
		}
	}
	value.push_back(infinity);
	// delta[(2*j + to)*mask_count + m]: the change of value[m] when entry j moves to haplotype "to"
	vector<unsigned int> delta(2 * n * mask_count);
	for (size_t j = 0; j < n; ++j) {
		const unsigned int shift = 2 * entry_individual[j];
		for (unsigned int to = 0; to < 2; ++to) {
			for (size_t m = 0; m < mask_count; ++m) {
				delta[(2 * j + to) * mask_count + m] = entry_cost[j][(masks[m] >> (shift + to)) & 1] - entry_cost[j][(masks[m] >> (shift + 1 - to)) & 1];
			}
		}
	}
	// assignment[i*width + a]: the masks allowed for configuration i, repeated up to the same width
	size_t width = 1;
	for (size_t i = 0; i < transmission_configurations; ++i) {
		width = std::max(width, assignments[i].size());
	}
	vector<unsigned int> assignment(transmission_configurations * width);
	for (size_t i = 0; i < transmission_configurations; ++i) {
		for (size_t a = 0; a < width; ++a) {
			assignment[i * width + a] = assignments[i].empty() ? mask_count : assignments[i][a % assignments[i].size()];
		}
	}

	// recombination[i*T + j]: the recombination cost between configuration j (previous column) and i
	vector<unsigned int> recombination(transmission_configurations * transmission_configurations);
	for (size_t i = 0; i < transmission_configurations; ++i) {
		for (size_t j = 0; j < transmission_configurations; ++j) {
			recombination[i * transmission_configurations + j] = popcount(i ^ j) * recombcost[column_index];
		}
	}
	// the change of the forward projection index when entry j changes its haplotype
	const vector<unsigned int>* forward_projection_mask = current_indexer->get_forward_projection_mask();
	vector<unsigned int> forward_flip(n, 0);
	if (forward_projection_mask != nullptr) {
		for (size_t j = 0; j < n; ++j) {
			if (forward_projection_mask->at(j) != (unsigned int)-1) {
				forward_flip[j] = ((unsigned int)1) << forward_projection_mask->at(j);
			}
		}
	}

	pedigree_column_t column;
	column.configurations = transmission_configurations;
	column.entries = n;
	column.mask_count = mask_count;
	column.width = width;
	column.backward_width = (column_index > 0) ? current_indexer->get_backward_projection_width() : 0;
	column.previous = (column_index > 0) ? &previous_projection_column->at(0, 0) : nullptr;
	column.projection = (current_projection_column != nullptr) ? &current_projection_column->at(0, 0) : nullptr;
	column.index_backtrace = (current_projection_column != nullptr) ? &index_backtrace_column->at(0, 0) : nullptr;
	column.transmission_backtrace = (current_projection_column != nullptr) ? &transmission_backtrace_column->at(0, 0) : nullptr;
	pedigree_optimum_t optimum = {optimal_score, optimal_score_index, optimal_transmission_value, previous_transmission_value};
	switch (transmission_configurations) {
	case 1:
		pedigree_column<1>(column, value, delta, assignment, recombination, forward_flip, optimum);
		break;
	case 4:
		pedigree_column<4>(column, value, delta, assignment, recombination, forward_flip, optimum);
		break;
	case 16:
		pedigree_column<16>(column, value, delta, assignment, recombination, forward_flip, optimum);
		break;
	default:
		pedigree_column<0>(column, value, delta, assignment, recombination, forward_flip, optimum);
	}
	optimal_score = optimum.score;
	optimal_score_index = optimum.index;
	optimal_transmission_value = optimum.transmission_value;
	previous_transmission_value = optimum.previous_transmission_value;
}


unsigned int PedigreeDPTable::get_optimal_score() {
	//if (backtrace_table.empty()) throw runtime_error("Empty backtrace table");
	return optimal_score;
}


void PedigreeDPTable::get_super_reads(std::vector<ReadSet*>* output_read_set, vector<unsigned int>* transmission_vector) {
	assert(output_read_set != nullptr);
	assert(output_read_set->size() == pedigree->size());
	assert(transmission_vector != nullptr);
	transmission_vector->clear();

	input_column_iterator.jump_to_column(0);
	const vector<unsigned int>* positions = input_column_iterator.get_positions();

	std::vector<std::pair<Read*,Read*>> superreads;
	for (unsigned int i=0; i<pedigree->size(); i++) {
		 superreads.emplace_back(
			new Read("superread_0_"+std::to_string(i), -1, -1, pedigree->index_to_id(i)),
			new Read("superread_1_"+std::to_string(i), -1, -1, pedigree->index_to_id(i))
		);
	}

	if (index_backtrace_table.empty()) {
		assert(!input_column_iterator.has_next());
	} else {
		// run through the file again with the input_column_iterator
		unsigned int i = 0; // column index
		while (input_column_iterator.has_next()) {
			const index_and_inheritance_t& v = index_path[i];
			unique_ptr<vector<const Entry *> > column = input_column_iterator.get_next();
			PedigreeColumnCostComputer cost_computer(*column, i, read_sources, pedigree, *pedigree_partitions[v.inheritance_value], distrust_genotypes);
			cost_computer.set_partitioning(v.index);

			auto population_alleles = cost_computer.get_alleles();
			
			// TODO: compute proper weights based on likelihoods.
			for (unsigned int k=0; k<pedigree->size(); k++) {
				superreads[k].first->addVariant(positions->at(i), population_alleles[k].allele0, population_alleles[k].quality);
				superreads[k].second->addVariant(positions->at(i), population_alleles[k].allele1, population_alleles[k].quality);
			}
			transmission_vector->push_back(v.inheritance_value);
			++i; // next column
		}
	}
	for(unsigned int k=0;k<pedigree->size();k++) {
		assert(output_read_set->at(k) != nullptr);
		output_read_set->at(k)->add(superreads[k].first);
		output_read_set->at(k)->add(superreads[k].second);
	}
}


vector<bool>* PedigreeDPTable::get_optimal_partitioning() {
	vector<bool>* partitioning = new vector<bool>(read_set->size(),false);

	for(size_t i=0; i< index_path.size(); ++i) {
		unsigned int mask = 1; // mask to pass over the partitioning (i.e., index)
		for(size_t j=0; j< indexers[i]->get_read_ids()->size(); ++j) {
			unsigned int index = index_path[i].index;
			if((index & mask) == 0) { // id at this index is in p0 (i.e., in the part.)
				partitioning->at(indexers[i]->get_read_ids()->at(j)) = true;
			}
			mask = mask << 1;
		}
	}
	
	return partitioning;
}
