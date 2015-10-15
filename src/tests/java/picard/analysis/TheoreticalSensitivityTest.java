/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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

import org.testng.Assert;
import org.testng.annotations.AfterTest;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Arrays;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.File;
import java.util.Scanner;

/**
 * Created by davidben on 5/18/15.
 */
public class TheoreticalSensitivityTest {

    private final static File TEST_DIR = new File("testdata/picard/analysis/TheoreticalSensitivity/");
    private final File depth = new File(TEST_DIR, "Solexa332667_DepthDist.histo");
    private final File baseQ = new File(TEST_DIR, "Solexa332667_BaseQ.histo");

    @Test
    public void testRouletteWheel() throws Exception {

        //test that a deterministic roulette wheel only gives one value
        final double[] deterministicWeights = {0.0, 1.0, 0.0};
        final TheoreticalSensitivity.RouletteWheel deterministicWheel = new TheoreticalSensitivity.RouletteWheel(deterministicWeights);
        for (int n = 0; n < 10; n++) Assert.assertEquals(deterministicWheel.draw(), 1);

        //test the sums of this deterministic wheel: a sum of n 1's equals n
        final List<ArrayList<Integer>> deterministicSums = deterministicWheel.sampleCumulativeSums(10, 1);
        for (int n = 0; n < 10; n++) Assert.assertEquals(deterministicSums.get(n).get(0), (Integer) n);

        //test that an unfair coin w/ probability 1/4 of heads gives roughly 1/4 heads
        /*final double p = 0.25;
        int total_heads = 0;
        final int N = 10000;
        final double stdDev = Math.sqrt(N*p*(1-p));   //the stddev of the sample total heads
        final List<Double> unfairCoinWeights = Arrays.asList(1-p, p);
        final TheoreticalSensitivity.RouletteWheel coinFlipWheel = new TheoreticalSensitivity.RouletteWheel(unfairCoinWeights);
        for (int n = 0; n < N; n++) total_heads += coinFlipWheel.draw();
        Assert.assertEquals(total_heads, N*p, 10*stdDev);
        */
    }

    @Test
    public void testProportionsAboveThresholds() throws Exception {
        final List<ArrayList<Integer>> sums = new ArrayList<ArrayList<Integer>>();
        sums.add(new ArrayList<Integer>(Arrays.asList(0,0,0)));
        sums.add(new ArrayList<Integer>(Arrays.asList(10, 10)));
        sums.add(new ArrayList<Integer>(Arrays.asList(5, 11, -2, 4)));
        final List<Double> thresholds = Arrays.asList(-1.0, 1.0, 6.0);
        Assert.assertEquals(sums.size(), 3);
        Assert.assertEquals(thresholds.size(), 3);

        final List<ArrayList<Double>> proportions = TheoreticalSensitivity.proportionsAboveThresholds(sums, thresholds);
        Assert.assertEquals(proportions.size(), 3);

        Assert.assertEquals(proportions.get(0).get(0), (double) 3/3);
        Assert.assertEquals(proportions.get(0).get(1), (double) 0/3);
        Assert.assertEquals(proportions.get(0).get(2), (double) 0/3);
        Assert.assertEquals(proportions.get(1).get(0), (double) 2/2);
        Assert.assertEquals(proportions.get(1).get(1), (double) 2/2);
        Assert.assertEquals(proportions.get(1).get(2), (double) 2/2);
        Assert.assertEquals(proportions.get(2).get(0), (double) 3/4);
        Assert.assertEquals(proportions.get(2).get(1), (double) 3/4);
        Assert.assertEquals(proportions.get(2).get(2), (double) 1/4);
    }

    @Test
    public void testHetAltDepthDistribution() throws Exception {
        final int N = 6;
        final double p = 0.5;
        final List<ArrayList<Double>> distribution = TheoreticalSensitivity.hetAltDepthDistribution(N);

        for (int n = 0; n < N-1; n++) {
            for (int m = 0; m <= n; m++) {
                //java has no built-in binomial coefficient
                //when this is in hellbender, use apache commons
                int binomialCoefficient = 1;
                for (int i = n; i > (n - m); i--) binomialCoefficient *= i;
                for (int i = m; i > 0; i--) binomialCoefficient /= i;

                Assert.assertEquals(distribution.get(n).get(m), binomialCoefficient*Math.pow(p,n));
            }
        }
    }

    //test that large-sample sums from the RouletteWheel converge to a normal distribution
    //using the empirical CDF as measured by proportionsAboveThresholds
    @Test
    public void testCentralLimitTheorem() throws Exception {
        //use a RouletteWheel that gives 0, 1, 2 with equal probability
        final double[] weights = {1.0, 1.0, 1.0};
        final TheoreticalSensitivity.RouletteWheel wheel = new TheoreticalSensitivity.RouletteWheel(weights);

        final int sampleSize = 1000;
        final int numSummands = 100;

        //the mean and standard deviation of a single roulette draw and of many draws
        final double muSingleDraw = 1.0;
        final double sigmaSingleDraw = Math.sqrt(2.0 / 3.0);
        final double mu = numSummands * muSingleDraw;
        final double sigma = Math.sqrt(numSummands) * sigmaSingleDraw;

        //test the sums of this deterministic wheel: a sum of n 1's equals n
        final List<ArrayList<Integer>> sums = wheel.sampleCumulativeSums(numSummands, sampleSize);
        //we only want the last set of sums, those with numSummands summands
        sums.subList(0, sums.size() - 1).clear();

        Assert.assertEquals(sums.size(), 1);

        //test whether the number of elements within one standard deviation agrees with the normal distribution
        final List<Double> thresholds = Arrays.asList(mu - sigma, mu + sigma);

        //sums is 1 x sampleSize, thresholds is a 2-vector, so proportions is 1 x 2
        final List<ArrayList<Double>> proportions = TheoreticalSensitivity.proportionsAboveThresholds(sums, thresholds);
        final double empiricalProportionWithinOneSigma = proportions.get(0).get(0) - proportions.get(0).get(1);

        //the proportion within one sigma for the normal distribution
        //hence whether any element falls within one sigma is a Bernoulli variable
        final double theoreticalProportionWithinOneSigma = 0.682689492;
        final double samplingStandardDeviationOfProportion = Math.sqrt(theoreticalProportionWithinOneSigma*(1-theoreticalProportionWithinOneSigma) /  sampleSize);

        Assert.assertEquals(empiricalProportionWithinOneSigma, theoreticalProportionWithinOneSigma, 5*samplingStandardDeviationOfProportion);
    }

    //Put it all together for deterministic quality and depths
    @Test
    public void testDeterministicQualityAndDepth() throws Exception {
        final double logOddsThreshold = 0.0;
        final double tolerance = 0.001;
        final int sampleSize = 1; //quality is deterministic, hence no sampling error
        for (int q = 5; q < 10; q++) {
            for (int n = 5; n < 10; n++) {
                final double minAltCount = 10*n*Math.log10(2)/q;  //alts required to call when log odds ratio threshold = 1
                double expectedResult = 0.0;

                final List<ArrayList<Double>> altCountProbabilities = TheoreticalSensitivity.hetAltDepthDistribution(n+1);
                for (int altCount = n; altCount > minAltCount; altCount--) {
                    expectedResult += altCountProbabilities.get(n).get(altCount);
                }

                //deterministic weights that always yield q are 0.0 for 0 through q - 1 and 1.0 for q
                final double[] qualityDistribution = new double[q+1];
                Arrays.fill(qualityDistribution, 0L);
                qualityDistribution[qualityDistribution.length-1]=1L;
                final double[] depthDistribution = new double[n+1];
                Arrays.fill(depthDistribution, 0L);
                depthDistribution[depthDistribution.length-1]=1L;

                final double result = TheoreticalSensitivity.hetSNPSensitivity(depthDistribution, qualityDistribution, sampleSize, logOddsThreshold);
                Assert.assertEquals(result, expectedResult, tolerance);
            }
        }
    }

    @Test
    public void testHetSensWGS() throws Exception {
        //Expect theoretical sens to be close to .9306 for Solexa-332667
        final double tolerance = 0.001;
        final double expectedResult = .9306;
        final int max = 250;
        final double [] depthDistribution = new double[max+1];
        final double [] qualityDistribution = new double[25];

        final Scanner scanDepth = new Scanner(depth);
        int i = 0;
        while (scanDepth.hasNextDouble()) {
            depthDistribution[i] = scanDepth.nextDouble();
            i++;
        }
        final Scanner scanBaseQ = new Scanner(baseQ);
        int j = 0;
        while (scanBaseQ.hasNextDouble()) {
            qualityDistribution[j] = scanBaseQ.nextDouble();
            j++;
        }

        final int sampleSize = 10000;
        final double logOddsThreshold = 3.0;

        final double result = TheoreticalSensitivity.hetSNPSensitivity(depthDistribution, qualityDistribution, sampleSize, logOddsThreshold);
        Assert.assertEquals(result, expectedResult, tolerance);
    }
}
