import pytest

from whatshap.core import Read, Variant, NumericSampleIds
from whatshap.variants import merge_two_reads, merge_reads
from whatshap.cli import PhasedInputReader
from whatshap.vcf import VcfReader


@pytest.mark.parametrize("merge", [merge_two_reads, merge_reads])
def test_merge_pair_without_shared_positions(merge):
    empty1 = Read("Name1")
    empty2 = Read("Name2")
    assert merge(empty1, empty2).name == "Name1"
    assert merge(empty2, empty1).name == "Name2"

    # add_variant parameters are: (position, allele, quality)
    left = Read("Name1")
    left.add_variant(100, 0, 31)
    left.add_variant(200, 0, 32)
    right = Read("Name2")
    right.add_variant(300, 1, 41)
    right.add_variant(400, 1, 42)

    expected = [
        Variant(100, 0, 31),
        Variant(200, 0, 32),
        Variant(300, 1, 41),
        Variant(400, 1, 42),
    ]
    assert expected == list(merge(left, right))
    assert expected == list(merge(right, left))

    outer = Read("Name1")
    outer.add_variant(100, 0, 31)
    outer.add_variant(400, 1, 42)
    inner = Read("Name2")
    inner.add_variant(200, 0, 32)
    inner.add_variant(300, 1, 41)
    assert expected == list(merge(inner, outer))
    assert expected == list(merge(outer, inner))


@pytest.mark.parametrize("merge", [merge_two_reads, merge_reads])
def test_merge_pair_with_shared_positions(merge):
    left = Read("Name1")
    left.add_variant(100, 0, 31)
    left.add_variant(200, 0, 32)
    left.add_variant(300, 0, 33)
    right = Read("Name2")
    right.add_variant(200, 0, 41)  # alleles disagree
    right.add_variant(300, 1, 42)  # alleles agree
    right.add_variant(400, 1, 43)

    expected = [
        Variant(100, 0, 31),
        Variant(200, 0, 32 + 41),
        Variant(300, 1, 42),
        Variant(400, 1, 43),
    ]
    assert expected == list(merge(left, right))
    assert expected == list(merge(right, left))


def test_merge_many_reads():
    reads = [
        Read("Name1"),
        Read("Name2"),
        Read("Name3"),
    ]
    reads[0].add_variant(100, 0, 31)
    reads[0].add_variant(200, 1, 32)
    reads[0].add_variant(300, 0, 33)

    reads[1].add_variant(200, 1, 41)
    reads[1].add_variant(400, 0, 42)
    reads[1].add_variant(500, 0, 43)

    reads[2].add_variant(200, 0, 51)
    reads[2].add_variant(500, 0, 52)
    reads[2].add_variant(600, 0, 53)

    expected = [
        Variant(100, 0, 31),
        Variant(200, 1, 73),  # see note: this depends on order of reads
        Variant(300, 0, 33),
        Variant(400, 0, 42),
        Variant(500, 0, 43 + 52),
        Variant(600, 0, 53),
    ]
    assert expected == list(merge_reads(*reads))

    # TODO merging should not depend on the order of reads
    expected[1] = Variant(200, 0, 51)
    assert expected == list(merge_reads(*reads[::-1]))


def test_allele_dection_01():
    path = "tests/data/alleledetection.biallelic"
    bam_reader = PhasedInputReader(
        [f"{path}.01.bam"],
        reference=None,
        numeric_sample_ids=NumericSampleIds(),
        ignore_read_groups=True,
        only_snvs=False,
        mapq_threshold=20,
    )
    vcf_reader = VcfReader(f"{path}.vcf", phases=False, only_snvs=False)
    sample = vcf_reader.samples[0]
    table = next(iter(vcf_reader))
    chromosome = table.chromosome
    readset, _vcf_source_ids = bam_reader.read(chromosome, table.variants, sample)
    expected = {}
    expected["Read01"] = [(102, 0), (105, 0)]
    expected["Read02"] = [(102, 0), (105, 1)]
    expected["Read03"] = [(102, 1), (105, 1)]
    expected["Read04"] = [(102, 0), (105, 1)]
    expected["Read05"] = [(102, 0), (105, 1)]
    expected["Read06"] = [(102, 0)]
    for read in readset:
        assert expected[read.name] == [(v.position, v.allele) for v in read]


def test_allele_dection_02():
    path = "tests/data/alleledetection.biallelic"
    bam_reader = PhasedInputReader(
        [f"{path}.02.bam"],
        reference=None,
        numeric_sample_ids=NumericSampleIds(),
        ignore_read_groups=True,
        only_snvs=False,
        mapq_threshold=20,
    )
    vcf_reader = VcfReader(f"{path}.vcf", phases=False, only_snvs=False)
    sample = vcf_reader.samples[0]
    table = next(iter(vcf_reader))
    chromosome = table.chromosome
    readset, _vcf_source_ids = bam_reader.read(chromosome, table.variants, sample)
    expected = {}
    expected["Read11"] = [(105, 0), (108, 0)]
    expected["Read12"] = [(105, 0), (108, 1)]
    expected["Read13"] = [(105, 0), (108, 0)]
    expected["Read14"] = [(105, 0), (108, 1)]
    expected["Read15"] = [(105, 0), (108, 1)]
    expected["Read16"] = [(105, 0), (108, 0)]
    expected["Read17"] = [(105, 0), (108, 0)]
    for read in readset:
        assert expected[read.name] == [(v.position, v.allele) for v in read]


def test_allele_dection_03():
    path = "tests/data/alleledetection.biallelic"
    bam_reader = PhasedInputReader(
        [f"{path}.03.bam"],
        reference=None,
        numeric_sample_ids=NumericSampleIds(),
        ignore_read_groups=True,
        only_snvs=False,
        mapq_threshold=20,
    )
    vcf_reader = VcfReader(f"{path}.vcf", phases=False, only_snvs=False)
    sample = vcf_reader.samples[0]
    table = next(iter(vcf_reader))
    chromosome = table.chromosome
    readset, _vcf_source_ids = bam_reader.read(chromosome, table.variants, sample)
    expected = {}
    expected["Read20"] = [(111, 0), (112, 0), (114, 0)]
    expected["Read21"] = [(111, 0), (112, 0), (114, 1)]
    expected["Read22"] = [(111, 1), (112, 0), (114, 0)]
    expected["Read23"] = [(111, 1), (112, 1), (114, 0)]
    expected["Read24"] = [(111, 0), (112, 0), (114, 0)]
    expected["Read25"] = [(111, 1), (112, 0), (114, 1)]
    expected["Read26"] = [(111, 1), (114, 1)]
    expected["Read27"] = [(117, 0)]
    expected["Read28"] = [(117, 1)]
    expected["Read29"] = []
    for read in readset:
        assert expected[read.name] == [(v.position, v.allele) for v in read]


def test_allele_dection_04():
    path = "tests/data/alleledetection.biallelic"
    bam_reader = PhasedInputReader(
        [f"{path}.04.bam"],
        reference=None,
        numeric_sample_ids=NumericSampleIds(),
        ignore_read_groups=True,
        only_snvs=False,
        mapq_threshold=20,
    )
    vcf_reader = VcfReader(f"{path}.vcf", phases=False, only_snvs=False)
    sample = vcf_reader.samples[0]
    table = next(iter(vcf_reader))
    chromosome = table.chromosome
    readset, _vcf_source_ids = bam_reader.read(chromosome, table.variants, sample)
    expected = {}
    expected["Read31"] = [(121, 0), (123, 0), (124, 0), (126, 0), (128, 0)]
    expected["Read32"] = [(121, 1), (123, 0), (124, 0), (126, 0), (128, 0)]
    expected["Read33"] = [(123, 0), (124, 0), (126, 0), (128, 0)]
    expected["Read34"] = [(121, 0), (123, 0), (124, 0), (126, 0), (128, 0)]
    expected["Read35"] = [(121, 0), (123, 0), (126, 0), (128, 0)]
    expected["Read36"] = [(121, 0), (123, 1), (124, 0), (126, 0), (128, 0)]
    expected["Read37"] = [(121, 0), (123, 1), (124, 0), (126, 1), (128, 0)]
    for read in readset:
        assert expected[read.name] == [(v.position, v.allele) for v in read]


def test_allele_dection_05():
    path = "tests/data/alleledetection.biallelic"
    for ref in [None, "tests/data/alleledetection.fasta"]:
        bam_reader = PhasedInputReader(
            [f"{path}.05.bam"],
            reference=ref,
            numeric_sample_ids=NumericSampleIds(),
            ignore_read_groups=True,
            only_snvs=False,
            mapq_threshold=20,
        )
        vcf_reader = VcfReader(f"{path}.vcf", phases=False, only_snvs=False)
        sample = vcf_reader.samples[0]
        table = next(iter(vcf_reader))
        chromosome = table.chromosome
        readset, _vcf_source_ids = bam_reader.read(chromosome, table.variants, sample)
        expected = {}
        expected["Read41"] = [(202, 0), (205, 0)]
        expected["Read42"] = [(202, 1), (205, 0)]
        expected["Read43"] = [(202, 0), (205, 1)]
        expected["Read44"] = [(202, 0), (205, 0)]
        expected["Read45"] = [(202, 0), (205, 1)]
        expected["Read46"] = [(202, 0)]
        expected["Read47"] = [(208, 0)]
        expected["Read48"] = [] if ref is None else [(208, 0)]
        expected["Read49"] = [] if ref is None else [(208, 0)]
        expected["Read50"] = [(208, 1)]
        for read in readset:
            assert expected[read.name] == [(v.position, v.allele) for v in read]


def test_allele_dection_multi_01():
    path = "tests/data/alleledetection.multiallelic"
    bam_reader = PhasedInputReader(
        [f"{path}.01.bam"],
        reference=None,
        numeric_sample_ids=NumericSampleIds(),
        ignore_read_groups=True,
        only_snvs=False,
        mapq_threshold=20,
    )
    vcf_reader = VcfReader(f"{path}.vcf", phases=False, only_snvs=False, mav=True)
    sample = vcf_reader.samples[0]
    table = next(iter(vcf_reader))
    chromosome = table.chromosome
    readset, _vcf_source_ids = bam_reader.read(chromosome, table.variants, sample)
    expected = {}
    expected["Read61"] = [(102, 0), (106, 0)]
    expected["Read62"] = [(102, 1), (106, 0)]
    expected["Read63"] = [(102, 1), (106, 2)]
    expected["Read64"] = [(102, 2), (106, 3)]
    for read in readset:
        assert expected[read.name] == [(v.position, v.allele) for v in read]


def test_allele_dection_multi_02():
    path = "tests/data/alleledetection.multiallelic"
    bam_reader = PhasedInputReader(
        [f"{path}.01.bam"],
        reference="tests/data/alleledetection.fasta",
        numeric_sample_ids=NumericSampleIds(),
        ignore_read_groups=True,
        only_snvs=False,
        mapq_threshold=20,
    )
    vcf_reader = VcfReader(f"{path}.vcf", phases=False, only_snvs=False, mav=True)
    sample = vcf_reader.samples[0]
    table = next(iter(vcf_reader))
    chromosome = table.chromosome
    readset, _vcf_source_ids = bam_reader.read(chromosome, table.variants, sample)
    expected = {}
    expected["Read61"] = [(102, 0), (106, 0)]
    expected["Read62"] = [(102, 1), (106, 0)]
    expected["Read63"] = [(102, 1), (106, 2)]
    expected["Read64"] = [(102, 2), (106, 3)]
    for read in readset:
        assert expected[read.name] == [(v.position, v.allele) for v in read]


def realign_reference(variants, restricted_genotypes, j, bam_read, reference, overhang):
    """The loop of ReadSetReader.detect_alleles_by_alignment over ReadSetReader.realign"""
    from whatshap._variants import _iterate_cigar
    from whatshap.variants import ReadSetReader

    cigartuples = bam_read.cigartuples
    result = []
    for index, i, consumed, query_pos in _iterate_cigar(variants, j, bam_read, cigartuples):
        allele, quality = ReadSetReader.realign(
            variants[index],
            restricted_genotypes[index] if restricted_genotypes else None,
            bam_read,
            cigartuples,
            i,
            consumed,
            query_pos,
            reference,
            overhang,
            False,
            None,
            None,
            None,
            False,
            None,
            7,
            40,
            25,
            None,
            None,
        )
        if allele in range(len(variants[index].get_alt_allele_list()) + 1):
            result.append((index, allele, quality))
    return result


def random_realignment_case(rng):
    import pysam
    from whatshap.vcf import BiallelicVcfVariant, MultiallelicVcfVariant

    bases = "ACGT"

    def seq(n):
        return "".join(rng.choice(bases) for _ in range(n))

    def other(b):
        return rng.choice([c for c in bases if c != b])

    reference = seq(300)
    # short repeats make indels ambiguous
    for _ in range(4):
        p = rng.randrange(20, 260)
        repeat = (seq(rng.choice([1, 2, 3])) * 10)[: rng.randrange(4, 12)]
        reference = reference[:p] + repeat + reference[p + len(repeat) :]
    variants = []
    pos = 15
    while True:
        pos += rng.randrange(1, 25)
        if pos >= 280:
            break
        kind = rng.random()
        if kind < 0.45:
            variants.append(BiallelicVcfVariant(pos, reference[pos], other(reference[pos])))
        elif kind < 0.65:
            variants.append(
                BiallelicVcfVariant(pos, reference[pos], reference[pos] + seq(rng.randrange(1, 6)))
            )
        elif kind < 0.85:
            variants.append(
                BiallelicVcfVariant(pos, reference[pos : pos + rng.randrange(2, 7)], reference[pos])
            )
        elif kind < 0.95:
            variants.append(
                MultiallelicVcfVariant(
                    pos, reference[pos], [other(reference[pos]), reference[pos] + seq(2)]
                )
            )
        else:
            variants.append(BiallelicVcfVariant(pos, reference[pos], "<DEL>"))
    by_position = {v.position: v for v in variants}

    start = rng.randrange(0, 60)
    end = rng.randrange(240, 300)
    ops = []  # (operation, reference length, read bases)
    if rng.random() < 0.3:
        ops.append(("S", 0, seq(rng.randrange(1, 15))))
    p = start
    while p < end:
        v = by_position.get(p)
        alt = v.get_alt_allele_list()[0] if v is not None else "<"
        if v is not None and not alt.startswith("<") and rng.random() < 0.5:
            r = v.reference_allele
            if len(r) == 1 and len(alt) == 1:
                ops.append(("X", 1, alt))
                p += 1
            elif len(alt) > len(r):
                ops.append(("=", 1, r))
                ops.append(("I", 0, alt[1:]))
                p += 1
            else:
                ops.append(("=", 1, r[0]))
                ops.append(("D", len(r) - 1, ""))
                p += len(r)
            continue
        e = rng.random()
        if e < 0.015:
            ops.append(("X", 1, other(reference[p])))
            p += 1
        elif e < 0.025:
            ops.append(("I", 0, seq(rng.randrange(1, 4))))
        elif e < 0.035 and p + 3 < end:
            n = rng.randrange(1, 3)
            ops.append(("D", n, ""))
            p += n
        else:
            ops.append(("=", 1, reference[p]))
            p += 1
    if rng.random() < 0.3:
        ops.append(("S", 0, seq(rng.randrange(1, 15))))
    use_m = rng.random() < 0.3
    cigar = []
    query = ""
    for op, n, b in ops:
        if use_m and op in "=X":
            op = "M"
        length = len(b) if op in "SIMX=" else n
        query += b
        if cigar and cigar[-1][0] == op:
            cigar[-1][1] += length
        else:
            cigar.append([op, length])
    bam_read = pysam.AlignedSegment()
    bam_read.query_name = "read"
    bam_read.query_sequence = query
    bam_read.reference_start = start
    bam_read.cigarstring = "".join(f"{n}{op}" for op, n in cigar)
    j = 0
    while j < len(variants) and variants[j].position < start:
        j += 1
    return reference, variants, bam_read, j


def test_detect_alleles_by_realignment_identical():
    # The default realignment gives the same alleles as the loop over ReadSetReader.realign
    import random
    from whatshap.core import Genotype
    from whatshap.variants import ReadSetReader

    rng = random.Random(1)
    observations = 0
    for case in range(300):
        reference, variants, bam_read, j = random_realignment_case(rng)
        restricted = None
        if case % 3 == 0:
            restricted = [Genotype([0, len(v.get_alt_allele_list())]) for v in variants]
        for overhang in (10, 3):
            expected = realign_reference(variants, restricted, j, bam_read, reference, overhang)
            detected = list(
                ReadSetReader.detect_alleles_by_alignment(
                    variants, restricted, j, bam_read, reference, overhang
                )
            )
            assert detected == expected
            observations += len(expected)
    assert observations > 10000
