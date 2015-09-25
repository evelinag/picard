/*
 * The MIT License
 *
 * Copyright (c) 2013 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.analysis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Metrics;
import picard.util.RExecutor;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Calculates and reports QC metrics for RRBS data based on the methylation status at individual C/G bases as well
 * as CpG sites across all reads in the input BAM/SAM file.
 *
 * @author jgentry@broadinstitute.org
 */

@CommandLineProgramProperties(
        usage = CollectRrbsMetrics.USAGE_SUMMARY + CollectRrbsMetrics.USAGE_DETAILS,
        usageShort = CollectRrbsMetrics.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class CollectRrbsMetrics extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Collects metrics about data generated with reduced representation bisulfite sequencing (RRBS).  ";
    static final String USAGE_DETAILS = "This tool calculates and reports QC metrics for RRBS data, based on the " +
            "methylation status of cytosine (C) bases in both CpG and non-CpG sites across all reads of a " +
            "BAM/SAM file.<br /><br />" +
            "" +
            "Cytosine methylation is a key component in epigenetic regulation of gene expression and frequently occurs at " +
            "CpG sites throughout the genome.  Bisulfite sequencing is a technique used to analyze the genome-wide " +
            "methylation profiles on a single nucleotide level [doi:10.1093/nar/gki901].   Sodium bisulfite efficiently and selectively deaminates " +
            "unmethylated cytosine residues to uracil without affecting 5-methyl cytosine (methylated).  " +
            "Using restriction enzymes and PCR to enrich for regions of the genome that have high CpG " +
            "content, the resulting reduced genome comprises ~1% of the original genome but includes key regulatory" +
            " sequences as well as repeated regions.<br /><br />" +
            "" +
            "The protocol involves several steps.  First, genomic DNA is digested with a restriction endonuclease such " +
            "as MspI, which targets CG dinucleotides.  This results in DNA fragments with CG at the ends.  Next, the " +
            "fragments are size selected (via gel electrophoresis), which facilitates the enrichment of CpG-containing sequences. " +
            " The next step involves bisulfite treatment to convert unmethylated C nucleotides to uracil (U), while methylated cytosines will remain intact." +
            "  The bisulfite-treated DNA is amplified with a proofreading-deficient DNA polymerase to facilitate amplification of both methylated cytosines as well as " +
            "the C -> U converted bases.  Subsequent to PCR amplification, each original unmethylated cytosine will be converted to either a T (+ strand) or an A (- strand), " +
            "while methylated C will remain a C (+ strand) or a G (- strand).  The PCR products are then sequenced" +
            " using conventional methods and aligned to a reference. <br /><br /> " +
            "" +
            "Since cytosine methylation is not exclusive for CpG \"hotspots\", the CollectRrbsMetrics tool outputs a summary table " +
            "indicating the number of CpG and non-CpG cytosines as well as their conversion C -> T (+ strand) or G -> A (- strand) rates. " +
            "  The tool also outputs the numbers of reads having no CpG sites, " +
            "and the numbers of reads discarded from the analysis due to inadequate size or excessive numbers of mismatches.<br /><br />" +

            "The tool also provides a table containing detailed information on CpG occurrence frequency, CpG conversion frequencies [C -> T (+ strand) or G -> A (- strand)], and the " +
            "specific locations of the CpG sites in the genome.   The conversion frequency helps determines the methylation status of a CpG site.<br /><br />" +
            "" +
            "Finally, the tool provides graphical representation of four metrics in the form of a \".pdf\" document." +
            "These metrics are the bisulfite conversion rate for CpG and non-CpG cytosines, a distribution of the numbers of CpG sites as a function of CpG conversion rate," +
            " the distribution of CpG sites by read coverage, and the numbers of reads discarded due to high numbers of mismatches or inadequate read size." +
            "" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar CollectRrbsMetrics \\<br />" +
            "      I=input.bam \\<br />" +
            "      M=rrbs_metrics \\<br />" +
            "      R=reference_sequence.fasta" +
            "</pre>" +
            "<hr />" +
            "" +
            "For additional details see: " +
            "http://broadinstitute.github.io/picard/picard-metric-definitions.html#RrbsCpgDetailMetrics and " +
            "http://broadinstitute.github.io/picard/picard-metric-definitions.html#RrbsSummaryMetrics";

// Path to R file for plotting purposes

private static final String R_SCRIPT = "picard/analysis/rrbsQc.R";

    @Option(doc = "The BAM or SAM file containing aligned reads. Must be coordinate sorted", shortName = StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;
    @Option(doc = "Base name for output files", shortName = StandardOptionDefinitions.METRICS_FILE_SHORT_NAME)
    public String METRICS_FILE_PREFIX;
    @Option(doc = "The reference sequence fasta file", shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME)
    public File REFERENCE;
    @Option(doc = "Minimum read length")
    public int MINIMUM_READ_LENGTH = 5;
    @Option(doc = "Threshold for base quality of a C base before it is considered")
    public int C_QUALITY_THRESHOLD = 20;
    @Option(doc = "Threshold for quality of a base next to a C before the C base is considered")
    public int NEXT_BASE_QUALITY_THRESHOLD = 10;
    @Option(doc = "Maximum percentage of mismatches in a read for it to be considered, with a range of 0-1")
    public double MAX_MISMATCH_RATE = 0.1;
    @Option(doc = "Set of sequence names to consider, if not specified all sequences will be used", optional = true)
    public Set<String> SEQUENCE_NAMES = new HashSet<String>();
    @Option(shortName = StandardOptionDefinitions.ASSUME_SORTED_SHORT_NAME,
            doc = "If true, assume that the input file is coordinate sorted even if the header says otherwise.")
    public boolean ASSUME_SORTED = false;
    @Option(shortName = "LEVEL", doc = "The level(s) at which to accumulate metrics.  ")
    public Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = CollectionUtil.makeSet(MetricAccumulationLevel.ALL_READS);

    public static final String DETAIL_FILE_EXTENSION = "rrbs_detail_metrics";
    public static final String SUMMARY_FILE_EXTENSION = "rrbs_summary_metrics";
    public static final String PDF_FILE_EXTENSION = "rrbs_qc.pdf";

    private static final Log log = Log.getInstance(CollectRrbsMetrics.class);

    public static void main(final String[] args) {
        new CollectRrbsMetrics().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        if (!METRICS_FILE_PREFIX.endsWith(".")) {
            METRICS_FILE_PREFIX = METRICS_FILE_PREFIX + ".";
        }
        final File SUMMARY_OUT = new File(METRICS_FILE_PREFIX + SUMMARY_FILE_EXTENSION);
        final File DETAILS_OUT = new File(METRICS_FILE_PREFIX + DETAIL_FILE_EXTENSION);
        final File PLOTS_OUT = new File(METRICS_FILE_PREFIX + PDF_FILE_EXTENSION);
        assertIoFiles(SUMMARY_OUT, DETAILS_OUT, PLOTS_OUT);

        final SamReader samReader = SamReaderFactory.makeDefault().open(INPUT);
        if (!ASSUME_SORTED && samReader.getFileHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            throw new PicardException("The input file " + INPUT.getAbsolutePath() + " does not appear to be coordinate sorted");
        }

        final ReferenceSequenceFileWalker refWalker = new ReferenceSequenceFileWalker(REFERENCE);
        final ProgressLogger progressLogger = new ProgressLogger(log);

        final RrbsMetricsCollector metricsCollector = new RrbsMetricsCollector(METRIC_ACCUMULATION_LEVEL, samReader.getFileHeader().getReadGroups(),
                C_QUALITY_THRESHOLD, NEXT_BASE_QUALITY_THRESHOLD, MINIMUM_READ_LENGTH, MAX_MISMATCH_RATE);

        for (final SAMRecord samRecord : samReader) {
            progressLogger.record(samRecord);
            if (!samRecord.getReadUnmappedFlag() && !isSequenceFiltered(samRecord.getReferenceName())) {
                final ReferenceSequence referenceSequence = refWalker.get(samRecord.getReferenceIndex());
                metricsCollector.acceptRecord(samRecord, referenceSequence);
            }
        }
        metricsCollector.finish();
        final MetricsFile<RrbsMetrics, Comparable<?>> rrbsMetrics = getMetricsFile();
        metricsCollector.addAllLevelsToFile(rrbsMetrics);

        // Using RrbsMetrics as a way to get both of the metrics objects through the MultiLevelCollector. Once
        // we get it out split it apart to the two separate MetricsFiles and write them to file
        final MetricsFile<RrbsSummaryMetrics, ?> summaryFile = getMetricsFile();
        final MetricsFile<RrbsCpgDetailMetrics, ?> detailsFile = getMetricsFile();
        for (final RrbsMetrics rrbsMetric : rrbsMetrics.getMetrics()) {
            summaryFile.addMetric(rrbsMetric.getSummaryMetrics());
            for (final RrbsCpgDetailMetrics detailMetric : rrbsMetric.getDetailMetrics()) {
                detailsFile.addMetric(detailMetric);
            }
        }
        summaryFile.write(SUMMARY_OUT);
        detailsFile.write(DETAILS_OUT);
        RExecutor.executeFromClasspath(R_SCRIPT, DETAILS_OUT.getAbsolutePath(), SUMMARY_OUT.getAbsolutePath(), PLOTS_OUT.getAbsolutePath());
        CloserUtil.close(samReader);
        return 0;
    }

    private boolean isSequenceFiltered(final String sequenceName) {
        return (SEQUENCE_NAMES != null) && (SEQUENCE_NAMES.size() > 0) && (!SEQUENCE_NAMES.contains(sequenceName));
    }

    private void assertIoFiles(final File summaryFile, final File detailsFile, final File plotsFile) {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(REFERENCE);
        IOUtil.assertFileIsWritable(summaryFile);
        IOUtil.assertFileIsWritable(detailsFile);
        IOUtil.assertFileIsWritable(plotsFile);
    }

    @Override
    protected String[] customCommandLineValidation() {
        final List<String> errorMsgs = new ArrayList<String>();
        if (MAX_MISMATCH_RATE < 0 || MAX_MISMATCH_RATE > 1) {
            errorMsgs.add("MAX_MISMATCH_RATE must be in the range of 0-1");
        }

        if (C_QUALITY_THRESHOLD < 0) {
            errorMsgs.add("C_QUALITY_THRESHOLD must be >= 0");
        }

        if (NEXT_BASE_QUALITY_THRESHOLD < 0) {
            errorMsgs.add("NEXT_BASE_QUALITY_THRESHOLD must be >= 0");
        }

        if (MINIMUM_READ_LENGTH <= 0) {
            errorMsgs.add("MINIMUM_READ_LENGTH must be > 0");
        }

        return errorMsgs.size() == 0 ? null : errorMsgs.toArray(new String[errorMsgs.size()]);
    }
}
