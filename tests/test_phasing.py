import random

from pytest import fixture
from whatshap.core import (
    Read,
    ReadSet,
    PedigreeDPTable,
    Pedigree,
    NumericSampleIds,
    PhredGenotypeLikelihoods,
    HapChatCore,
)
from whatshap.testhelpers import (
    string_to_readset,
    brute_force_phase,
    canonic_index_to_biallelic_gt,
    canonic_index_list_to_biallelic_gt_list,
)


@fixture(params=["whatshap", "hapchat"])
def algorithm(request):
    return request.param


def test_phase_empty_readset(algorithm):
    rs = ReadSet()
    recombcost = [1, 1]
    genotypes = canonic_index_list_to_biallelic_gt_list([1, 1])
    pedigree = Pedigree(NumericSampleIds())
    genotype_likelihoods = [None, None]
    pedigree.add_individual("individual0", genotypes, genotype_likelihoods)

    if algorithm == "hapchat":
        dp_table = HapChatCore(rs)
    else:
        dp_table = PedigreeDPTable(rs, recombcost, pedigree)

    _ = dp_table.get_super_reads()


def compare_phasing_brute_force(
    superreads, cost, partition, readset, all_heterozygous, weights=None, algorithm="whatshap"
):
    """Compares DPTable based phasing to brute force phasing and returns string representation of superreads."""
    assert len(superreads) == 2
    assert len(superreads[0]) == len(superreads[1])
    for v1, v2 in zip(*superreads):
        assert v1.position == v2.position
    haplotypes = tuple(sorted("".join(str(v.allele) for v in sr) for sr in superreads))
    (
        expected_cost,
        expected_partition,
        solution_count,
        expected_haplotype1,
        expected_haplotype2,
    ) = brute_force_phase(readset, all_heterozygous)
    print()
    print(superreads[0])
    print(superreads[1])
    print("Partition:", partition)
    print("Expected: ", expected_partition)
    print("Haplotypes:")
    print(haplotypes[0])
    print(haplotypes[1])
    print("Expected Haplotypes:")
    print(expected_haplotype1)
    print(expected_haplotype2)
    print("Cost:", cost)
    print("Expected cost:", expected_cost)
    # TODO: implement the reporting of an optimal partitioning in hapchat
    if algorithm == "whatshap":
        inverse_partition = [1 - p for p in partition]
        assert (partition == expected_partition) or (inverse_partition == expected_partition)
    assert solution_count == 1
    assert cost == expected_cost
    assert (haplotypes == (expected_haplotype1, expected_haplotype2)) or (
        haplotypes == (expected_haplotype2, expected_haplotype1)
    )


def check_phasing_single_individual(reads, algorithm="whatshap", weights=None):
    # 0) set up read set
    readset = string_to_readset(reads, weights)
    positions = readset.get_positions()

    # for hapchat
    if algorithm == "hapchat":
        dp_table = HapChatCore(readset)
        superreads = dp_table.get_super_reads()
        cost = dp_table.get_optimal_cost()
        partition = dp_table.get_optimal_partitioning()
        compare_phasing_brute_force(
            superreads[0][0], cost, partition, readset, True, weights, algorithm
        )
        return

    # 1) Phase using PedMEC code for single individual
    for all_heterozygous in [False, True]:
        recombcost = [1] * len(positions)  # recombination costs 1, should not occur
        pedigree = Pedigree(NumericSampleIds())
        genotype_likelihoods = [
            None if all_heterozygous else PhredGenotypeLikelihoods([0, 0, 0])
        ] * len(positions)
        pedigree.add_individual(
            "individual0",
            [canonic_index_to_biallelic_gt(1) for i in range(len(positions))],
            genotype_likelihoods,
        )  # all genotypes heterozygous
        dp_table = PedigreeDPTable(
            readset, recombcost, pedigree, distrust_genotypes=not all_heterozygous
        )
        superreads, transmission_vector = dp_table.get_super_reads()
        cost = dp_table.get_optimal_cost()
        # TODO: transmission vectors not returned properly, see issue 73
        assert len(set(transmission_vector)) == 1
        partition = dp_table.get_optimal_partitioning()
        compare_phasing_brute_force(
            superreads[0], cost, partition, readset, all_heterozygous, weights
        )

    # 2) Phase using PedMEC code for trios with two "empty" individuals (i.e. having no reads)
    for all_heterozygous in [False, True]:
        recombcost = [1] * len(positions)  # recombination costs 1, should not occur
        pedigree = Pedigree(NumericSampleIds())
        genotype_likelihoods = [
            None if all_heterozygous else PhredGenotypeLikelihoods([0, 0, 0])
        ] * len(positions)
        pedigree.add_individual(
            "individual0",
            [canonic_index_to_biallelic_gt(1) for _ in range(len(positions))],
            genotype_likelihoods,
        )  # all genotypes heterozygous
        pedigree.add_individual(
            "individual1",
            [canonic_index_to_biallelic_gt(1) for _ in range(len(positions))],
            genotype_likelihoods,
        )  # all genotypes heterozygous
        pedigree.add_individual(
            "individual2",
            [canonic_index_to_biallelic_gt(1) for _ in range(len(positions))],
            genotype_likelihoods,
        )  # all genotypes heterozygous
        pedigree.add_relationship("individual0", "individual1", "individual2")
        dp_table = PedigreeDPTable(
            readset, recombcost, pedigree, distrust_genotypes=not all_heterozygous
        )
        cost = dp_table.get_optimal_cost()
        superreads, transmission_vector = dp_table.get_super_reads()
        assert len(set(transmission_vector)) == 1
        partition = dp_table.get_optimal_partitioning()
        compare_phasing_brute_force(
            superreads[0], cost, partition, readset, all_heterozygous, weights
        )


def test_phase_trivial(algorithm):
    reads = """
          11
           01
        """
    check_phasing_single_individual(reads, algorithm)


def test_phase1(algorithm):
    reads = """
     10
     010
     010
    """
    check_phasing_single_individual(reads, algorithm)


def test_phase2(algorithm):
    reads = """
      1  11010
      00 00101
      001 0101
    """
    check_phasing_single_individual(reads, algorithm)


def test_phase3(algorithm):
    reads = """
      1  11010
      00 00101
      001 01010
    """
    check_phasing_single_individual(reads, algorithm)


def test_phase4(algorithm):
    reads = """
      1  11010
      00 00101
      001 01110
       1    111
    """
    check_phasing_single_individual(reads, algorithm)


def test_phase4a(algorithm):
    reads = """
      1  11010
      00 00101
      001 01110
       1    111
    """
    check_phasing_single_individual(reads, algorithm)


# note: these final two tests do not apply to hapchat because their
# (brute force solutions) are weighted phasings -- hapchat does not do
# weights (for the time being)
def test_phase5():
    reads = """
      0             0
      110111111111
      00100
           0001000000
           000
            10100
                  101
    """
    check_phasing_single_individual(reads)


def test_weighted_phasing1():
    reads = """
      1  11010
      00 00101
      001 01110
       1    111
    """
    weights = """
      2  13112
      11 23359
      223 56789
       2    111
    """
    check_phasing_single_individual(reads, "whatshap", weights)


def random_readset(rng, n_positions, n_reads, max_span, qualities, max_first=None):
    positions = [10 * (i + 1) for i in range(n_positions)]
    reads = []
    for index in range(n_reads):
        first = rng.randrange(max_first or (n_positions - 1))
        last = min(n_positions - 1, first + rng.randrange(1, max_span))
        read = Read(f"read{index}", 50, 0)
        covered = [i for i in range(first, last + 1) if i in (first, last) or rng.random() < 0.8]
        for i in covered:
            read.add_variant(positions[i], rng.randrange(2), rng.choice(qualities))
        reads.append(read)
    readset = ReadSet()
    for read in reads:
        readset.add(read)
    readset.sort()
    return readset


def phase_single_individual(readset, allow_single_individual_dp):
    positions = readset.get_positions()
    pedigree = Pedigree(NumericSampleIds())
    pedigree.add_individual(
        "individual0",
        [canonic_index_to_biallelic_gt(1) for _ in positions],
        [None] * len(positions),
    )
    dp_table = PedigreeDPTable(
        readset,
        [1] * len(positions),
        pedigree,
        positions=positions,
        allow_single_individual_dp=allow_single_individual_dp,
    )
    superreads, transmission_vector = dp_table.get_super_reads()
    haplotypes = [[(v.position, v.allele, v.quality) for v in sr] for sr in superreads[0]]
    return (
        haplotypes,
        transmission_vector,
        dp_table.get_optimal_cost(),
        dp_table.get_optimal_partitioning(),
    )


def test_single_individual_dp_identical():
    # The DP for one heterozygous individual gives the same result as the general DP,
    # also for ties (few distinct qualities) and for columns of more than 16 reads
    rng = random.Random(1)
    for n_positions, n_reads, max_span, qualities, max_first, repeats in [
        (2, 3, 2, [1], None, 20),
        (30, 40, 8, [1, 2], None, 20),
        (60, 60, 12, [1, 2, 3, 30], None, 20),
        (200, 120, 20, [5, 10, 20, 30, 40, 60], None, 20),
        (12, 18, 12, [1, 2, 3], 2, 3),  # columns of more than 16 active reads
    ]:
        for _ in range(repeats):
            readset = random_readset(rng, n_positions, n_reads, max_span, qualities, max_first)
            assert phase_single_individual(readset, True) == phase_single_individual(
                readset, False
            )
