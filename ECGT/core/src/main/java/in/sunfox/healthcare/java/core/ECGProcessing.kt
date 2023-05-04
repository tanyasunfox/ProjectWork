package `in`.sunfox.healthcare.java.core

import kotlin.math.abs
import kotlin.math.ceil
import kotlin.math.roundToInt
import kotlin.math.sqrt



    class ECGProcessing(
        private var ecgPoints: ArrayList<Double>,
        private var leadPosition: Int ,
        private var applyFilter: Boolean,
        private var adjustRPeaks: Boolean,
        private var ecgPointsForAugmentedLead: ArrayList<Double>? = null
    ) {


        lateinit var prStartIndices: ArrayList<Double>
        var tWaveEndPoints = ArrayList<Double>()
            private set
        var qrsStopIndices = ArrayList<Double>()
            private set

        var sPoints = ArrayList<Double>()
            private set

        var qPoints = ArrayList<Double>()
            private set

        var tPoints = ArrayList<Double>()
            private set

        var qTcInterval = 0.0 // added For New Algorithm: change is in return type Integer to Double
            get() {
                if (!isSignalProcessed) return 0.0
                field = if (heartRate in 60..100) {
//                (DecimalFormat("#.##").format(qtInterval.toDouble() / sqrt(rrInterval.toDouble() / 1000.0))).toDouble()

                    (Filters.localDecimalFormat((qtInterval.toDouble() / sqrt(rrInterval.toDouble() / 1000.0)),"#.##")).toDouble()
                } else {
//                DecimalFormat("#.##").format(qtInterval + 1000.0 * (0.154 * (1 - rrInterval.toDouble() / 1000.0))).toDouble()
                    Filters.localDecimalFormat((qtInterval + 1000.0 * (0.154 * (1 - rrInterval.toDouble() / 1000.0))),"#.##").toDouble()
                }
                return field
            }

        var heartRate = 0
            get() {
                return if (!isSignalProcessed) 0
                else field
            }

        var msDifference: Double = if(ecgPoints.size<1000000)
        {
            if(ecgPoints.size>3000)
                10000.0 / (ecgPoints.size-8)
            else
                5000.0 / (ecgPoints.size-8)
        }
        else
            300000.0 / (ecgPoints.size-8) //1.98

        var rPeakPoints: ArrayList<Double> = ArrayList()
            private set
        var rrIntervals: ArrayList<Double> = ArrayList()
            private set
        var qRSIntervals: ArrayList<Double> = ArrayList()
            private set
        private var qtIntervals: ArrayList<Double> = ArrayList()
        var prIntervals: ArrayList<Double> = ArrayList()
        var pPoints: ArrayList<Double> = ArrayList()
            private set
        var prStopIndices: ArrayList<Double> = ArrayList()
            private set
        private var isSignalProcessed = false

        var sampleToIgnore=arrayListOf<Int>()

        var stSegmentInMv = 0.0
        var averagePAmplitudeInLead = 0.0 //added for new Algorithm
        var averageQAmplitudeInLead = 0.0 //added for new Algorithm
        var averageSAmplitudeInLead = 0.0 //added for new Algorithm
        var averageTAmplitudeInLead = 0.0 //added for new Algorithm
        var averageRAmplitudeInLead = 0.0 //added for new Algorithm
        var pWidth = 0.0 //added for new Algorithm
        var tWidth = 0.0 //added for new Algorithm
        var qrsDirectionUpward = false //added for new Algorithm
        var valueRS = 0.0 //added for new Algorithm
        var ventricularActivationLOR = 0.0 //added for new Algorithm
        var ventricularActivationROR = 0.0 //added for new Algorithm
        var concavity = false //added for new Algorithm
        var frequencyOfPatternInQRS=0  //added for new Algorithm
        var frequencyOfPatternInRR=0 //added for new Algorithm
        var pAmplitudeArrayInMv = arrayListOf<Double>()  //added for new Algorithm
        var TRRatioSatisfy =false
        var TSRatioSatisfy=false

        /////////////////////////////////////////////////////////
        var tentedTWave = false
        var multiPWave=false
        var pPeaksIndex = arrayListOf<Double>()
        var prIntervalInMilliSec = arrayListOf<Double>()
        var arrayRateOfMultiPWave = arrayListOf<Double>()
        var pWavePresent=true
        /////////////////////////////////////////////////////////


        private var unFilteredEcgPoints = arrayListOf<Double>()

        private var isAugmentedLead: Boolean = false

        fun process(
            rPointsForAug: ArrayList<Double>? = null,
            sPointsForAug: ArrayList<Double>? = null,
            qPointsForAug: ArrayList<Double>? = null,
            pPointsForAug: ArrayList<Double>? = null,
            tPointsForAug: ArrayList<Double>? = null
        ) {
            isSignalProcessed = true
            for(i in ecgPoints.indices){
                sampleToIgnore.add(0);
            }
            if (rPointsForAug == null) {
                rPeakPoints = calculateRPoints(ecgPoints)
                if (adjustRPeaks)
                    rPeakPoints = adjustRPeaks(unFilteredEcgPoints, rPeakPoints)
            } else {
                rPeakPoints = rPointsForAug
            }

            calculateRRInterval()
            calculateHeartRate()
            processForPQST(
                sPointsForAug,
                qPointsForAug,
                pPointsForAug,
                tPointsForAug
            )
        }

        init {
            isAugmentedLead =
                leadPosition >= 8 && leadPosition <= 11

            ecgPoints.forEach {
                unFilteredEcgPoints.add(it)
            }
            if (applyFilter)
                ecgPoints = Filters.movingAverage(ecgPoints)
            // msDifference = 10000.0 / ecgPoints.size
        }

        fun adjustRPeaks(
            signalOnConsideration: ArrayList<Double>,
            seudoPoint: ArrayList<Double>
        ): ArrayList<Double> {
            val pointRSample: ArrayList<Double> = ArrayList()
            var i = 0
            var checkRightLeft = 0
            while (i < seudoPoint.size) {
                var j = seudoPoint[i]
                var maximum = signalOnConsideration[j.toInt()]
                if (signalOnConsideration[seudoPoint[i].toInt()] - signalOnConsideration[seudoPoint[i].toInt() - 1] > 0 && signalOnConsideration[seudoPoint[i].toInt()] - signalOnConsideration[seudoPoint[i].toInt() + 1] > 0) {
                    pointRSample.add(seudoPoint[i]);
                } else {

                    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    if (signalOnConsideration[seudoPoint[i].toInt()] - signalOnConsideration[seudoPoint[i].toInt() - 1] > 0) {
                        checkRightLeft = 1
                        if (seudoPoint[i] + 50 <= signalOnConsideration.size) {
                            checkRightLeft *= 50
                        } else checkRightLeft *= (signalOnConsideration.size - seudoPoint[i].toInt() - 1)
                    } else {
                        checkRightLeft = -1
                        checkRightLeft = if (seudoPoint[i] - 50 >= 1) checkRightLeft * 50 else 0
                    }

                    var maximumSample = j
                    if (checkRightLeft > 0) {
                        while (j < seudoPoint[i] + checkRightLeft) {
                            if (signalOnConsideration[j.toInt()] > signalOnConsideration[j.toInt() - 1] && signalOnConsideration[j.toInt()] > signalOnConsideration[j.toInt() + 1] || signalOnConsideration[j.toInt()] > signalOnConsideration[j.toInt() - 1] && signalOnConsideration[j.toInt()].equals(
                                    signalOnConsideration[j.toInt() + 1]
                                ) || signalOnConsideration[j.toInt()].equals(signalOnConsideration[j.toInt() - 1]) && signalOnConsideration[j.toInt()] > signalOnConsideration[j.toInt() + 1]
                            ) {
                                if (maximum < signalOnConsideration[j.toInt()]) {
                                    maximum = signalOnConsideration[j.toInt()]
                                    maximumSample = j
                                    j = seudoPoint[i] + checkRightLeft  // saturday
                                }
                            }
                            j += 1
                        }
                    } else {
                        while (j > seudoPoint[i] + checkRightLeft) {
                            if (signalOnConsideration[j.toInt()] > signalOnConsideration[j.toInt() - 1] && signalOnConsideration[j.toInt()] > signalOnConsideration[j.toInt() + 1] || signalOnConsideration[j.toInt()] > signalOnConsideration[j.toInt() - 1] && signalOnConsideration[j.toInt()].equals(
                                    signalOnConsideration[j.toInt() + 1]
                                ) || signalOnConsideration[j.toInt()].equals(signalOnConsideration[j.toInt() - 1]) && signalOnConsideration[j.toInt()] > signalOnConsideration[j.toInt() + 1]
                            ) {
                                if (maximum < signalOnConsideration[j.toInt()]) {
                                    maximum = signalOnConsideration[j.toInt()]
                                    maximumSample = j
                                    j = seudoPoint[i] + checkRightLeft  // saturday
                                }
                            }
                            j -= 1
                        }
                    }
                    pointRSample.add(maximumSample)
                    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
                }
                i += 1
            }
            return pointRSample
        }


        private fun processForPQST(
            sPointsForAug: ArrayList<Double>?,
            qPointsForAug: ArrayList<Double>?,
            pPointsForAug: ArrayList<Double>?,
            tPointsForAug: ArrayList<Double>?
        ) {

//        if (isAugmentedLead && applyFilter) {
//            this.ecgPoints = unFilteredEcgPoints;
//        }
            sPoints = sPointsForAug ?: calculateSPoints()
//        if(leadPosition==3)
//            println()

            for(i in 0 until sPoints.size) {
                if(ecgPoints.get(sPoints.get(i).toInt()) < 10.0) {
                    for(j in sPoints.get(i).toInt() + 1 until sPoints.get(i).toInt()+50) {
                        if(ecgPoints.get(j) >= 10.0) {
                            sPoints.set(i,j.toDouble())
                            break
                        }
                    }
                }

                //////////////////////////code to find s in pvc/////////////////////////////
                if(sPoints[i].toInt()==sampleToIgnore[sPoints[i].toInt()] /* || ecgPoints[sPoints[i].toInt()]<10**/) {
                    var control=10.0;
                    for(j in rPeakPoints[i].toInt() until sPoints[i].toInt()) {
                        if (((ecgPoints[rPeakPoints[i].toInt()] - ecgPoints[j]) / (j - rPeakPoints[i])) > ((ecgPoints[j] - ecgPoints[sPoints[i].toInt()]) / (sPoints[i] - j))) {
                            if (ecgPoints[j - 1] - ecgPoints[j] < control) {
                                sPoints[i] = j.toDouble();
                                control = ecgPoints[j - 1] - ecgPoints[j];
                            } else {
                                if (control != 10.0) {
                                    break;
                                }
                            }
                        }
                    }
                }
                //////////////////////////code to find s in pvc/////////////////////////////
            }

//        ////////////////Removing RS Pair not according to the definition//////////////////////
//        var holdArrayR=ArrayList<Int>();
//        var holdArrayS=ArrayList<Int>();
//        var count=0;
//        for(i in 0 until sPoints.size) {
//            if((ecgPoints[rPeakPoints[i].toInt()] - ecgPoints[sPoints[i].toInt()]) / ecgPoints[rPeakPoints[i].toInt()] < .1){
//            holdArrayR.add(i-count);
//            holdArrayS.add(i-count);
//            count+=1;
//
//            }
//        }
//
//        for(i in 0 until holdArrayR.size)
//            rPeakPoints.removeAt(holdArrayR.get(i))
//
//        for(i in 0 until holdArrayS.size)
//            sPoints.removeAt(holdArrayS.get(i))
//        ////////////////Removing RS Pair not according to the definition//////////////////////

            qPoints = qPointsForAug ?: calculateQpoints()
            pPoints = pPointsForAug ?: calculatePPoints(qPoints)
//        tPoints = tPointsForAug ?: calculateTpoints()
//        qPoints = calculateQpoints()
//        pPoints = calculatePPoints(qPoints)
//        tPoints = calculateTpoints(sPoints)
            prStartIndices = calculatePRstartIndices(pPoints)
            prStopIndices = calculatePRstopIndices(qPoints)
            qrsStopIndices = calculateQRSstopPoints(sPoints, prStopIndices)
            tPoints = tPointsForAug ?: calculateTpoints()
            tWaveEndPoints = calculateTwaveEndPoints(tPoints)
            if (leadPosition >= 0 && leadPosition <= 7) {
                stSegmentInMv = calculateSTSegmentInMv(qrsStopIndices, prStopIndices)
                averagePAmplitudeInLead = average(amplitudeArrayInMv(pPoints, prStopIndices, 1)); //added for new Algorithm
                pAmplitudeArrayInMv = amplitudeArrayInMv(pPoints, prStopIndices, 1); //added for new Algorithm
                averageQAmplitudeInLead = average(amplitudeArrayInMv(qPoints, prStopIndices, 1)); //added for new Algorithm
                averageSAmplitudeInLead = average(amplitudeArrayInMv(sPoints, prStopIndices, 0)); //added for new Algorithm
                averageTAmplitudeInLead = average(amplitudeArrayInMv(tPoints, prStopIndices, 0)); //added for new Algorithm
                averageRAmplitudeInLead =
                    average(amplitudeArrayInMv(rPeakPoints, prStopIndices, 0)); //added for new Algorithm
                pWidth = pWidth(pPoints, prStopIndices, ecgPoints) //added for new Algorithm
                tWidth = tWidth(tPoints, prStopIndices, ecgPoints) //added for new Algorithm
                qrsDirectionUpward = qrsDirectionUpward(rPeakPoints, sPoints, ecgPoints) //added for new Algorithm
                valueRS = ratioRS(rPeakPoints, sPoints, prStopIndices, "NonAugmented"); //added for new Algorithm
                concavity = concavity(sPoints, tPoints, prStopIndices, ecgPoints) //added for new Algorithm

            } else {
                stSegmentInMv = calculateAugmentedSTElevation(
                    qrsStopIndices,
                    prStopIndices,
                    ecgPointsForAugmentedLead!!
                )
                averagePAmplitudeInLead =
                    average(amplitudeArrayInMvForAug(pPoints, prStopIndices, 1)); //added for new Algorithm
                averageQAmplitudeInLead =
                    average(amplitudeArrayInMvForAug(qPoints, prStopIndices, 1)); //added for new Algorithm
                averageSAmplitudeInLead =
                    average(amplitudeArrayInMvForAug(sPoints, prStopIndices, 0)); //added for new Algorithm
                averageTAmplitudeInLead =
                    average(amplitudeArrayInMvForAug(tPoints, prStopIndices, 0)); //added for new Algorithm
                averageRAmplitudeInLead =
                    average(amplitudeArrayInMvForAug(rPeakPoints, prStopIndices, 0)); //added for new Algorithm
                pWidth = pWidth(pPoints, prStopIndices, ecgPointsForAugmentedLead!!) //added for new Algorithm
                tWidth = tWidth(tPoints, prStopIndices, ecgPointsForAugmentedLead!!) //added for new Algorithm
                qrsDirectionUpward =
                    qrsDirectionUpward(rPeakPoints, sPoints, ecgPointsForAugmentedLead!!) //added for new Algorithm
                valueRS = ratioRS(rPeakPoints, sPoints, prStopIndices, "Augmented"); //added for new Algorithm
                concavity =
                    concavity(sPoints, tPoints, prStopIndices, ecgPointsForAugmentedLead!!) //added for new Algorithm
            }
            prIntervals = calculatePRIntervals(prStartIndices, prStopIndices)
            qRSIntervals = calculateQRSIntervals(qrsStopIndices, prStopIndices)
            qtIntervals = calculateQTIntervals(tWaveEndPoints, prStopIndices)
            ventricularActivationLOR = ventricularActivationLOR(rPeakPoints, prStopIndices) //added for new Algorithm
            ventricularActivationROR = ventricularActivationROR(rPeakPoints, sPoints) //added for new Algorithm
            frequencyOfPatternInQRS = patternCheckInQRSArray()
            frequencyOfPatternInRR= patternCheckInRRIntervalArray();
            TRRatioSatisfy=tRAmp();
            TSRatioSatisfy=tSAmp();
            tentedTWavePresent();
            multiPwave();

//        if(leadPosition==7)
//            println()
        }


        val prInterval: Int
            get() {
                if (!isSignalProcessed) return 0
                var total = 0.0
                for (PR in prIntervals) total += PR
                return if (prIntervals.size > 0)
                    (total / prIntervals.size * msDifference).roundToInt()
                else 0;
            }

        val qrsInterval: Int
            get() {
                if (!isSignalProcessed) return 0
                var total = 0.0
                for (QRS in qRSIntervals) total += QRS
                return if (qRSIntervals.size > 0)
                    (total / qRSIntervals.size * msDifference).roundToInt()
                else 0;
            }

        val qtInterval: Int
            get() {
                if (!isSignalProcessed) return 0
                var total = 0.0
                for (QT in qtIntervals) total += QT
                return if (qtIntervals.size > 0)
                    (total / qtIntervals.size * msDifference).roundToInt()
                else 0;
            }


        val isEcgSignalCompatibleForProcessing: Boolean
            get() {
                val derivedSignal = calculateDerivative(ecgPoints)
                val maximaPoints = calculateMaximaPoints(derivedSignal)
                var isTheSignalECG=false;
//            var isTheSignalNotInRange=false;
//            for (i in 0 until ecgPoints.size)
//            {
//                if(ecgPoints[i]<995 || ecgPoints[i]>1000)
//                {
//                    isTheSignalNotInRange=true;
//                    break
//                }
//            }

                var high=0.0;
                var low=1000.0;
                for (i in 0 until ecgPoints.size)
                {
                    if(ecgPoints[i]>high)
                    {
                        high=ecgPoints[i];
                    }
                    if(ecgPoints[i]<low)
                    {
                        low=ecgPoints[i];
                    }
                }

//            if(maximaPoints.size > 1) {
                if (high - low > 30)
                    isTheSignalECG = true
//            }

                return isTheSignalECG
            }


        private fun calculateQTIntervals(
            tEnd: ArrayList<Double>,
            prStop: ArrayList<Double>?
        ): ArrayList<Double> {
            val qt = ArrayList<Double>()
            for (i in 1 until tEnd.size) qt.add(tEnd[i] - prStop!![i - 1])
            return qt
        }

        private fun calculateQRSIntervals(
            qrsStop: ArrayList<Double>,
            prStop: ArrayList<Double>?
        ): ArrayList<Double> {
            val qrs = ArrayList<Double>()
            for (i in 0 until prStop!!.size - 1) qrs.add(qrsStop[i] - 15 - prStop[i])  //added for new Algorithm: change is 30 to 15
            return qrs
        }

        public fun calculatePRIntervals(
            prStart: ArrayList<Double>,
            PRend: ArrayList<Double>?
        ): ArrayList<Double> {
            val prInts = ArrayList<Double>()
            for (i in prStart.indices) {
                prInts.add(PRend!![i] - prStart[i])
            }
            return prInts
        }

        private fun calculateSTSegmentInMv(
            qrsStop: ArrayList<Double>,
            prStopIndices: ArrayList<Double>
        ): Double {
            var totalSTE = 0.0
            if(prStopIndices.size!=0) {
                for (i in 0 until qrsStop.size) {
                    totalSTE += ecgPoints[qrsStop[i].toInt()] - ecgPoints[prStopIndices[i].toInt()]
                }
                if (qrsStop.size != 0)
                    return totalSTE * 0.0029296875 / qrsStop.size.toDouble()
                else
                    return 0.0
            }
            else
                return 0.0

        }

        fun calculateAugmentedSTElevation(
            qrsStop: ArrayList<Double>,
            prStopIndices: ArrayList<Double>,
            augmentedSignalData: ArrayList<Double>
        ): Double {
            var inversionFactor = 0;
            var totalSTE = 0.0
            if(prStopIndices.size!=0) {
                for (i in 0 until qrsStop.size) {
                    totalSTE += augmentedSignalData[qrsStop[i].toInt()] - augmentedSignalData[prStopIndices[i].toInt()]
//            println("${augmentedSignalData[qrsStop[i].toInt()]} &%& ${augmentedSignalData[prStopIndices!![i].toInt()]}")
                }
                var isInvertedSTElevation = false

                if (leadPosition == 10) {
                    for (i in 0 until rPeakPoints.size - 1) {
                        if (ecgPointsForAugmentedLead?.get(rPeakPoints[i].toInt())!! < ecgPointsForAugmentedLead?.get(
                                sPoints[i].toInt()
                            )!!
                        ) {
                            isInvertedSTElevation = true
//                   break
                        }
                    }
                }
                if (qrsStop.size != 0)
                    return totalSTE * 0.0029296875 * (if (isInvertedSTElevation) -1 else 1) / qrsStop.size.toDouble()
                else
                    return 0.0
            }
            else
                return 0.0
        }

        //================================================added for new Algorithm========================================================//

        private fun patternCheckInQRSArray(): Int
        {
            var isTherePattern=0
            for(i in 0 until qRSIntervals.size-1)
            {
                if(qRSIntervals.get(i)*.625>qRSIntervals.get(i+1) || qRSIntervals.get(i+1)*.625>qRSIntervals.get(i))
                    isTherePattern=isTherePattern+1;
            }
            return isTherePattern
        }
        private fun patternCheckInRRIntervalArray(): Int
        {
            var isTherePattern=0
            for(i in 0 until rrIntervals.size-1)
            {
                if(rrIntervals.get(i)>(1.75*rrIntervals.get(i+1)) || rrIntervals.get(i+1)>(1.75*rrIntervals.get(i)))
                    isTherePattern=isTherePattern+1;
            }
            return isTherePattern
        }
        private fun amplitudeArrayInMv(
            array: ArrayList<Double>,
            prStopIndices: ArrayList<Double>,
            whichParameter: Int
        ): ArrayList<Double> {
            val FACTOR = 0.0029296875
            val amplitudeArrayInMv = ArrayList<Double>();

            if (prStopIndices.size != 0) {
                if (array.size != 0) {
                    if (whichParameter == 1) {
                        for (i in 0 until array.size) {
                            amplitudeArrayInMv.add(
                                ((ecgPoints[array.get(i).toInt()] - ecgPoints[prStopIndices[i].toInt()]) * FACTOR)
                            )
                        }
                    } else {
                        for (i in 0 until array.size) {
                            if (i == 0)
                                amplitudeArrayInMv.add(
                                    ((ecgPoints[array.get(i).toInt()] - ecgPoints[prStopIndices[0].toInt()]) * FACTOR)
                                )
                            else
                                amplitudeArrayInMv.add(
                                    ((ecgPoints[array.get(i).toInt()] - ecgPoints[prStopIndices[i - 1].toInt()]) * FACTOR)
                                )
                        }
                    }
                }
            }

            return amplitudeArrayInMv
        }

        private fun amplitudeArrayInMvForAug(
            array: ArrayList<Double>,
            prStopIndices: ArrayList<Double>,
            whichParameter: Int
        ): ArrayList<Double> {
            val FACTOR = 0.0029296875
            var amplitudeArrayInMv = ArrayList<Double>();
            if (prStopIndices.size != 0) {
                if (array.size != 0) {
                    if (whichParameter == 1) {
                        for (i in 0 until array.size) {
                            amplitudeArrayInMv.add(
                                ((ecgPoints[array.get(i).toInt()] - ecgPoints[prStopIndices[i].toInt()]) * FACTOR)
                            )
                        }
                    } else {
                        for (i in 0 until array.size) {
                            if (i == 0)
                                amplitudeArrayInMv.add(
                                    ((ecgPointsForAugmentedLead?.get(
                                        array.get(i).toInt()
                                    )!! - (ecgPointsForAugmentedLead?.get(prStopIndices[0].toInt()))!!) * FACTOR)
                                )
                            else
                                amplitudeArrayInMv.add(
                                    ((ecgPointsForAugmentedLead?.get(
                                        array.get(i).toInt()
                                    )!! - (ecgPointsForAugmentedLead?.get(prStopIndices[i - 1].toInt()))!!) * FACTOR)
                                )
                        }
                    }
                }
            }

            return amplitudeArrayInMv
        }

        private fun pWidth(pPoint: ArrayList<Double>, prStop: ArrayList<Double>, signal: ArrayList<Double>): Double {
            val width = ArrayList<Double>()
            val toleranceFactor = .10
            var j: Int
            for (i in 0 until pPoint.size) {
                j = pPoint[i].toInt();
                while (j != -1) {
                    if (((signal[j] - signal[prStop[i].toInt()]) / (signal[pPoint[i].toInt()] - signal[prStop[i].toInt()])) < toleranceFactor
                    ) {
                        width.add(2.0 * (j.toDouble() - pPoint[i]));
                        j = -2;
                    }

                    if (j < signal.size - 1)
                        j += 1;
                    else
                        j = -1;
                }
            }

            return average(width) / 500.0
        }

        private fun tWidth(tPoint: ArrayList<Double>, prStop: ArrayList<Double>, signal: ArrayList<Double>): Double {
            var j: Int
            var first: Int
            var last: Int
            var width = 0.0

            for (i in 1 until tPoint.size) {
                first = 0;
                last = 0;
                j = tPoint[i].toInt();
                while (j != -1) {
                    if (((signal[j]) < (signal[prStop[i - 1].toInt()]) && (signal[tPoint[i].toInt()]) > (signal[prStop[i - 1].toInt()])) ||
                        ((signal[j]) > (signal[prStop[i - 1].toInt()]) && (signal[tPoint[i].toInt()]) < (signal[prStop[i - 1].toInt()]))
                    ) {
                        first = j;
                        j = 0;
                    }
                    if (j - 1 > 0)
                        j -= 1;
                    else
                        j = -1;
                }

                j = tPoint[i].toInt();
                while (j != -1) {
                    if (((signal[j]) < (signal[prStop[i - 1].toInt()]) && (signal[tPoint[i].toInt()]) > (signal[prStop[i - 1].toInt()])) ||
                        ((signal[j]) > (signal[prStop[i - 1].toInt()]) && (signal[tPoint[i].toInt()]) < (signal[prStop[i - 1].toInt()]))
                    ) {
                        last = j;
                        j = -2;
                    }
                    if (j + 1 < signal.size)
                        j += 1;
                    else
                        j = -1;
                }

                width += (last - first);
            }
            if(tPoint.size>1)
                return width / ((tPoint.size - 1) * 500.0)
            else
                return 0.0

        }

        private fun qrsDirectionUpward(
            rPeak: ArrayList<Double>,
            sPeak: ArrayList<Double>,
            signal: ArrayList<Double>
        ): Boolean {
            var count = 0

            for (i in 0 until rPeak.size - 1) {
                if (signal.get(rPeak.get(i).toInt()) > signal.get(sPeak.get(i).toInt())) count += 1
            }

            return count > (rPeak.size * .05)

        }

        private fun ratioRS(
            rPeak: ArrayList<Double>,
            sPeak: ArrayList<Double>,
            prStop: ArrayList<Double>,
            string: String
        ): Double {
            val valueRS = ArrayList<Double>()
            val rAmplitude: ArrayList<Double>
            val sAmplitude: ArrayList<Double>

            if (rPeak.size * sPeak.size != 0) {
                if (string == "NonAugmented") {
                    rAmplitude = amplitudeArrayInMv(rPeak, prStop, 0)
                    sAmplitude = amplitudeArrayInMv(sPeak, prStop, 0)
                } else {
                    rAmplitude = amplitudeArrayInMvForAug(rPeak, prStop, 0)
                    sAmplitude = amplitudeArrayInMvForAug(sPeak, prStop, 0)
                }

                for (i in 0 until rAmplitude.size - 1) {
                    if (sAmplitude[i] != 0.0) {
                        valueRS.add(rAmplitude[i] / sAmplitude[i])
                    }
                }
            }
            return average(valueRS)
        }

        private fun ventricularActivationLOR(rPeak: ArrayList<Double>, prStop: ArrayList<Double>): Double {
            var ventricularActivation = 0.0
            if (prStop.size != 0) {
                if (rPeak.size != 0) {
                    for (i in 1 until rPeak.size)
                        ventricularActivation += ((rPeak.get(i) - prStop.get(i - 1)) / 500.0);
                }
            }
            if (rPeak.size > 1)
                return ventricularActivation / (rPeak.size - 1).toDouble()
            else
                return 0.0
        }

        private fun ventricularActivationROR(rPeak: ArrayList<Double>, sPeak: ArrayList<Double>): Double {
            var ventricularActivation = 0.0
            if (rPeak.size * sPeak.size != 0) {
                for (i in 0 until rPeak.size - 1)
                    ventricularActivation += ((sPeak.get(i) - rPeak.get(i)) / 500.0);
            }
            var x = ventricularActivation / (rPeak.size - 1).toDouble();
            if (rPeak.size > 1)
                return ventricularActivation / (rPeak.size - 1).toDouble()
            else
                return 0.0
        }

        private fun concavity(
            sPeak: ArrayList<Double>,
            tPeak: ArrayList<Double>,
            prStop: ArrayList<Double>,
            signal: ArrayList<Double>
        ): Boolean {
            val width = 15.0
            var count: Double
            var filterSignal = ArrayList<Double>()
            val slope = ArrayList<Double>()
            var countConcave = 0
            var countConvex = 0
            var indicate = 0
            var q: Double

            if (sPeak.size * tPeak.size != 0) {
                if(sPeak[0]>2.5*width && (signal.size-tPeak.get(tPeak.size - 1).toInt()-1)>2.5*width) {
                    for (i in 0 until signal.size) {
                        if (i >= (sPeak.get(0).toInt()) - (width.toInt()) && i <= (tPeak.get(tPeak.size - 1)
                                .toInt()) + (width.toInt())
                        ) {
                            count = 0.0;
                            for (j in i - (width.toInt()) until i + (width.toInt()) + 1) {
                                count = count + signal.get(j);
                            }
                            filterSignal.add(count / (2 * width + 1));
                        } else
                            filterSignal.add(signal.get(i));
                    }
                }
                else
                    filterSignal=signal;

                var j = 0
                var last = 0
                for (i in 1 until sPeak.size-1) {
                    j = tPeak.get(i).toInt();
                    last = 0;
                    while (j != -1 && j >= sPeak.get(i).toInt()) {
                        if ((((signal.get(j) - signal.get(qrsStopIndices.get(i - 1).toInt())) <= .70 * (signal.get(
                                tPeak.get(i).toInt()
                            ) - signal.get(qrsStopIndices.get(i - 1).toInt()))) && 0 < signal[tPeak.get(i)
                                .toInt()] - signal[qrsStopIndices.get(i).toInt()])

                            || ((signal.get(j) - signal.get(qrsStopIndices.get(i - 1).toInt())) >= .70 * (signal.get(
                                tPeak.get(i).toInt()
                            ) - signal.get(qrsStopIndices.get(i - 1).toInt())) && 0 > signal[tPeak.get(i)
                                .toInt()] - signal[qrsStopIndices.get(i).toInt()])
                        ) {
                            last = j;
                            j = 0;
                        }
                        j -= 1;
                    }


                    slope.clear()
                    if (last != 0) {
                        q = (tPeak[i] - sPeak[i] + 1.0) / 3.0;
                        q = ceil(q);

                        j = qrsStopIndices[i].toInt();
                        while (j <= last) {
                            if (j + q.toInt() <= last)
                                slope.add((filterSignal[j + q.toInt()] - filterSignal[j]) / q);
                            else
                                slope.add((filterSignal[last] - filterSignal[j]) / (last - j));

                            if (j + q <= last)
                                j += q.toInt();
                            else
                                j = last + 1;
                        }

                        countConcave = 0
                        countConvex = 0

                        for (j in 1 until slope.size) {
                            if (((slope[j] - slope[j - 1]) > 0.0 && 0 < signal[tPeak.get(i)
                                    .toInt()] - signal[qrsStopIndices.get(i).toInt()])
                                || ((slope[j] - slope[j - 1]) < 0.0 && 0 > signal[tPeak.get(i)
                                    .toInt()] - signal[qrsStopIndices.get(i).toInt()])
                            )
                                countConcave += 1;
                            else
                                countConvex += 1;
                        }
                        if (countConcave >= countConvex)
                            indicate += 1;
                    }

                }
                if(sPeak.size>1)
                    return ((indicate.toDouble() / (sPeak.size - 1).toDouble()) > .5)
                else
                    return false
            } else
                return false
        }



        private fun average(array: ArrayList<Double>): Double {
            var average = 0.0
            for (i in 0 until array.size) {
                average += array[i];
            }
            return if (array.size != 0)
                average / array.size.toDouble()
            else
                0.0
        }

        private fun tRAmp():Boolean
        {
            var countHowManyQRSSatisfy=0
            var array=ArrayList<Double>()
            for(i in 1 until rPeakPoints.size-1)
            {
                array.add((ecgPoints[tPoints[i].toInt()]-ecgPoints[prStopIndices[i-1].toInt()])/(ecgPoints[rPeakPoints[i].toInt()]-ecgPoints[prStopIndices[i-1].toInt()]));
                if ((ecgPoints[tPoints[i].toInt()]-ecgPoints[prStopIndices[i-1].toInt()])/(ecgPoints[rPeakPoints[i].toInt()]-ecgPoints[prStopIndices[i-1].toInt()])<=-.3) {
                    countHowManyQRSSatisfy++
                }
            }
            return countHowManyQRSSatisfy>=2
        }

        private fun tSAmp():Boolean
        {
            var countHowManyQRSSatisfy=0
            for(i in 1 until rPeakPoints.size-1)
            {
                if ((ecgPoints[tPoints[i].toInt()]-ecgPoints[prStopIndices[i-1].toInt()])/(ecgPoints[sPoints[i].toInt()]-ecgPoints[prStopIndices[i-1].toInt()])<=-.3 && (ecgPoints[rPeakPoints[i].toInt()]-ecgPoints[prStopIndices[i].toInt()])*0.0029296875<.2) {
                    countHowManyQRSSatisfy++
                }
            }
            return countHowManyQRSSatisfy>=2
        }

        private fun multiPwave() {

            val middleOfRRPeaks = ArrayList<Double>()
            val lengthOfHalfRRPeaks = ArrayList<Double>()
            var i: Int
            var j: Int
            var q: Int
            var first: Int
            var check: Int
            var sign: Int
            var width: Int
            var countLeft: Int
            var countRight: Int
            var check1: Int
            var howManyZerosAreThere: Int
            var p: Int
            var temp: Boolean
            var indexOfP: Double
            var indexOfPRStart: Double
            var indexOfPRStop: Double
//        var zone = 0.0
            var countMultiPWavePair=0.0
            i = 1
            while (i < rPeakPoints.size) {
                middleOfRRPeaks.add(Math.ceil((rPeakPoints.get(i) + rPeakPoints.get(i - 1)) / 2.0))
                lengthOfHalfRRPeaks.add(Math.ceil((rPeakPoints.get(i) - rPeakPoints.get(i - 1)) / 2.0))
                i++
            }
            val arraySegment = java.util.ArrayList<Double>()
            val point = java.util.ArrayList<Double>()
            val multiPWave = java.util.ArrayList<Double>()
            i = 0
            while (i < middleOfRRPeaks.size) {
                arraySegment.clear()
                j = middleOfRRPeaks[i].toInt()
                while (j <= qPoints.get(i).toInt()) {
                    arraySegment.add(ecgPoints.get(j))
                    j++
                }
                point.clear()
                multiPWave.clear()
                point.add(qPoints.get(i))
                first = 1
                j = qPoints.get(i).toInt()
                check = 0
                sign = 1
                indexOfPRStop = 0.0
                while (j >= (middleOfRRPeaks[i] + sPoints.get(i)) / 2.0) {
                    width = if (check == 0) 30 else {
                        if (sign == 1) 20 else 10
                    }
                    var slopeAtPoint = 0.0
                    for (a in 1..10) {
                        slopeAtPoint = (ecgPoints.get(j) - ecgPoints.get(j - a)) / a
                        if (slopeAtPoint != 0.0) break
                    }
                    if (first == 1 && slopeAtPoint != 0.0) sign = (slopeAtPoint / Math.abs(slopeAtPoint)).toInt()
                    if (sign * ecgPoints.get(j) > sign * ecgPoints.get(j - 1) && sign * ecgPoints.get(j) > sign * ecgPoints.get(j + 1) || sign * ecgPoints.get(
                            j
                        ) == sign * ecgPoints.get(j - 1) && sign * ecgPoints.get(j) > sign * ecgPoints.get(j + 1) || sign * ecgPoints.get(
                            j
                        ) > sign * ecgPoints.get(j - 1) && sign * ecgPoints.get(j) == sign * ecgPoints.get(j + 1)
                    ) {
                        if (sign == 1 && indexOfPRStop == 0.0) indexOfPRStop = j.toDouble()
                        q = 1
                        countLeft = 0
                        countRight = 0
                        temp = false
                        if (sign == -1 && (qPoints.get(i) - j) / (qPoints.get(i) - middleOfRRPeaks[i]) < .1) temp =
                            true
                        howManyZerosAreThere = 0
                        while (q <= width && !temp) {
                            if (sign * ecgPoints.get(j) > sign * ecgPoints.get(j + q) || sign * ecgPoints.get(j) == sign * ecgPoints.get(
                                    j + q
                                )
                            ) countRight = countRight + 1
                            if (sign * ecgPoints.get(j) > sign * ecgPoints.get(j - q) || sign * ecgPoints.get(j) == sign * ecgPoints.get(
                                    j - q
                                )
                            ) {
                                countLeft = countLeft + 1
                                if (sign * ecgPoints.get(j) == sign * ecgPoints.get(j - q)) howManyZerosAreThere =
                                    howManyZerosAreThere + 1
                            }
                            q = q + 1
                        }
                        if (countLeft == width && countRight == width || countLeft == width && countRight > .8 * width || countRight == width && countLeft > .8 * width) {
                            point.add(j.toDouble())
                            j = j - howManyZerosAreThere
                            first = 0
                            sign = sign * -1
                            check = check + 1
                        }
                    }
                    j = j - 1
                }
                if (point.size > 2) {
                    point.add(Math.ceil((middleOfRRPeaks[i] + sPoints.get(i)) / 2))
                    indexOfP = 0.0
                    indexOfPRStart = 0.0
                    j = 1
                    while (j < point.size - 1) {
                        if (ecgPoints.get(point[j].toInt()) > ecgPoints.get(point[j - 1].toInt()) && ecgPoints.get(point[j].toInt()) > ecgPoints.get(
                                point[j + 1].toInt()
                            )
                        ) {
                            val ref1 = Math.ceil((qPoints.get(i) * 2 + point[j]) / 3)
                            val ref2 = 2 * point[j] - ref1
                            val ref: Double =
                                ecgPoints.get(ref2.toInt()) + (ecgPoints.get(ref1.toInt()) - ecgPoints.get(ref2.toInt())) * (point[j] - ref2) / (ref1 - ref2)
                            val value = (ecgPoints.get(point[j].toInt()) - ref) * 0.0029296875 * 10
                            if ((ecgPoints.get(point[j].toInt()) - ref) * 0.0029296875 * 10 > .4) {
                                if ((point[j] - middleOfRRPeaks[i]) / lengthOfHalfRRPeaks[i] > .17) {
                                    indexOfP = point[j]
                                    indexOfPRStart = 2.0 * point[j] - point[j - 1]
                                    p = point[j].toInt()
                                    while (p != -1) {
                                        if((ecgPoints.get(p) < ecgPoints.get(p - 1) && ecgPoints.get(p) < ecgPoints.get(p + 1)) ||
                                            (ecgPoints.get(p).equals(ecgPoints.get(p - 1)) && ecgPoints.get(p)<(ecgPoints.get(p + 1))) ||
                                            (ecgPoints.get(p) < ecgPoints.get(p - 1) && ecgPoints.get(p).equals(ecgPoints.get(p + 1)))){
                                            indexOfPRStart = p.toDouble()
                                            p = 0
                                        }
                                        p = p - 1
                                    }
                                    if (2 * point[j] - point[j - 1] > indexOfPRStart) indexOfPRStart =
                                        1.5 * point[j] - .5 * point[j - 1]
                                    if (indexOfPRStop == point[j]) indexOfPRStop = qPoints.get(i)
                                    multiPWave.add(point[j])
                                    check1 = 0
                                    var compare = 0.0
                                    q = j
                                    while (q < point.size - 2) {
                                        var refMultiP: Double
                                        if (Math.abs((ecgPoints.get(point[q + 1].toInt()) - ecgPoints.get(point[q - 1].toInt())) / (point[q + 1] - point[q - 1])) < 1 / Math.sqrt(
                                                3.0
                                            ) || q == j
                                        ) refMultiP =
                                            ecgPoints.get(
                                                point[q - 1].toInt()
                                            ) + (ecgPoints.get(point[q + 1].toInt()) - ecgPoints.get(point[q - 1].toInt())) * (point[q] - point[q - 1]) / (point[q + 1] - point[q - 1]) else refMultiP =
                                            Math.max(
                                                ecgPoints.get(
                                                    point[q - 1].toInt()
                                                ), ecgPoints.get(point[q + 1].toInt())
                                            )
                                        if (q == j) compare = ecgPoints.get(point[q].toInt()) - refMultiP
                                        if (q > j && ecgPoints.get(point[q].toInt()) - refMultiP >= compare * .9) {
                                            check1 = if (check1 == 0 && point[q] < middleOfRRPeaks[i]) break else {
                                                multiPWave.add(point[q])
                                                check1 + 1
                                            }
                                        }
                                        q += 2
                                    }
                                    if (multiPWave.size >= 2) {
                                        countMultiPWavePair++;
                                        var rate = 0.0
                                        for (qq in 0 until multiPWave.size - 1) rate =
                                            rate + multiPWave[qq] - multiPWave[qq + 1]
                                        arrayRateOfMultiPWave.add(60000 / (rate / (multiPWave.size - 1) * (10000.0 / ecgPoints.size)))
                                    }
                                    break
                                }
                            }
                        }
                        j++
                    }

                    if (indexOfP != 0.0) {
                        pPeaksIndex.add(indexOfP)
                        prIntervalInMilliSec.add((indexOfPRStop - indexOfPRStart) * (10000.0 / ecgPoints.size))
                    } else {
//                    zone = rPeakPoints.get(i)
                        pPeaksIndex.add(-999.0)
                    }
                } else {
//                zone = rPeakPoints.get(i)
                    pPeaksIndex.add(-999.0)
                }
                i = i + 1
            }
            var countThePeaks = 0
            i = 0
            while (i < pPeaksIndex.size) {
                if (pPeaksIndex.get(i) != -999.0) countThePeaks = countThePeaks + 1
                i++
            }

//        if(zone>3000.0) {
//            if (pPeaksIndex.size() - countThePeaks > 1) {
//                isPwavePresent = false;
//                isMultiPWave = false;
//            }
//        }
//           else {
//            if(pPeaksIndex.size() - countThePeaks > 0) {
//                isPwavePresent = false;
//                isMultiPWave = false;
//            }
//        }

//        int a=0;
//        int b=0;
//        if(zone>3000.0)
//        {
//            b=rPeakPoints.get(rPeakPoints.size()-1).intValue();
//            a=((int)zone)-(3000-b);
//        }
//        else {
//            b=2999;
//            a=0;
//        }
//        if(countThePeaks != pPeaksIndex.size())
            if(countMultiPWavePair>2)
                this.multiPWave = true

            if (pPeaksIndex.size - countThePeaks > 3) {
                pWavePresent = false
                this.multiPWave = false
            }
//        if(leadPosition==7)
//            println()

        }

        private fun tentedTWavePresent(){
            if(rPeakPoints.size!=0 && tPoints.size!=0 && prStopIndices.size!=0 && qrsStopIndices.size!=0) {

                var tAmplitude = amplitudeArrayInMv(tPoints, prStopIndices, 0);
                var rAmplitude = amplitudeArrayInMv(rPeakPoints, prStopIndices, 0);

                var right70 = ArrayList<Int>();
                var left70 = ArrayList<Int>();

                var referencePRSTOP = ArrayList<Double>();
                var referenceQRSTOP = ArrayList<Double>();

                if((rPeakPoints.get(1)-prStopIndices.get(0))<rPeakPoints.get(0))
                    referencePRSTOP.add(rPeakPoints.get(0)-(rPeakPoints.get(1)-prStopIndices.get(0)));
                else
                    referencePRSTOP.add(0.0);

                for(i in 0 until prStopIndices.size)
                {
                    referencePRSTOP.add(prStopIndices.get(i))
                }
                referenceQRSTOP.add(rPeakPoints.get(0)+(qrsStopIndices.get(0)-rPeakPoints.get(1)));
                for(i in 0 until qrsStopIndices.size)
                {
                    referenceQRSTOP.add(qrsStopIndices.get(i))
                }


                for(i in 0 until tPoints.size) {

                    var start = 0;
                    var j = tPoints.get(i).toInt();
                    while (j != -1) {
                        if (((tAmplitude.get(i) > 0 && ecgPoints.get(j) < ecgPoints.get(referencePRSTOP.get(i).toInt())) ||
                                    (tAmplitude.get(i) < 0 && ecgPoints.get(j) > ecgPoints.get(referencePRSTOP.get(i).toInt()))) &&
                            j >= referenceQRSTOP.get(i).toInt())
                        {
                            start = j + 1;
                            j = 0;
                        }
                        else
                        {
                            start = referenceQRSTOP.get(i).toInt();
                            j = 0;
                        }
                        j = j - 1;
                    }

//               var finish = 0;
//               j = tPoints.get(i).toInt();
//               while (j != -1) {
//                   if (((tAmplitude.get(i) > 0 && ecgPoints.get(j) < ecgPoints.get(referencePRSTOP.get(i).toInt())) ||
//                        (tAmplitude.get(i) < 0 && ecgPoints.get(j) > ecgPoints.get(referencePRSTOP.get(i).toInt()))) &&
//                         j <= 2 * tPoints.get(i).toInt() - referenceQRSTOP.get(i).toInt())
//                   {
//                       finish = j - 1;
//                       j = -2;
//                   }
//                   else
//                   {
//                       finish = 2 * tPoints.get(i).toInt() - referenceQRSTOP.get(i).toInt();
//                       j = -2;
//                   }
//                   j = j + 1;
//               }

                    j = tPoints.get(i).toInt();
                    while(j!=-1)
                    {
                        if((ecgPoints.get(j) - ecgPoints.get(prStopIndices.get(i).toInt())) / (ecgPoints.get(tPoints.get(i).toInt()) - ecgPoints.get(prStopIndices.get(i).toInt())) < .70 && j + 1 >= start) {
                            right70.add(j + 1);
                            j = 0;
                        }
                        else {
                            right70.add(start);
                            j = 0;
                        }
                        j = j - 1;
                    }

                    j = tPoints.get(i).toInt();
                    while(j!= -1) {
                        if((ecgPoints.get(j) - ecgPoints.get(prStopIndices.get(i).toInt()) / (ecgPoints.get(tPoints.get(i).toInt()) - ecgPoints.get(prStopIndices.get(i).toInt()))) < .70 && j - 1 <= 2 * tPoints.get(i) - start) {
                            if(j<=ecgPoints.size)
                                left70.add(j - 1);
                            else
                                left70.add(ecgPoints.size-1);

                            j = -2;
                        }
                        else {
                            if( (2 * tPoints.get(i).toInt() - start) < ecgPoints.size)
                                left70.add(2 * tPoints.get(i).toInt() - start);
                            else
                                left70.add(ecgPoints.size-1);

                            j = -2;
                        }
                        j = j + 1;
                    }

                }

                var countPeakness=0;
                var countSymmetricity=0;

                for(i in 0 until tPoints.size) {
                    if(right70.get(i)>0&&left70.get(i)>0) {
                        if (tAmplitude.get(i) >= rAmplitude.get(i) && tAmplitude.get(i) * 10 > 13) {
                            var width = left70.get(i) - right70.get(i);
                            var height = ecgPoints.get(tPoints.get(i).toInt()) - (ecgPoints.get(right70.get(i)) + ecgPoints.get(left70.get(i))) / 2;
                            if ((2 * height * 0.0029296875 * 10 * 500 / (width * 25)) > sqrt(3.0))
                                countPeakness = countPeakness + 1;

                            var width1 = tPoints.get(i).toInt() - right70.get(i);
                            var width2 = left70.get(i) - tPoints.get(i).toInt();

                            val maximumValue = if (width1 > width2) width1 else width2

                            if (abs((width1 - width2).toDouble()) / maximumValue <= .2)
                                countSymmetricity = countSymmetricity + 1

                        }
                    }

                }
                if(countPeakness.toDouble()/tPoints.size.toDouble()>=.42 && countSymmetricity.toDouble()/tPoints.size.toDouble()>=.42 )
                    tentedTWave=true;
//                    if(leadPosition==1)
//                        println()

            }
//tentedTWave=false;

        }

        //===============================================================================================================================//

        private fun calculateSPoints(): ArrayList<Double> {
            val sPoints = ArrayList<Double>()
            for (i in 0 until rPeakPoints.size - 1) {
                var j = rPeakPoints[i].toInt() + 10


                while(j>=rPeakPoints[i].toInt()+10) {
                    if(ecgPoints[j] < ecgPoints[rPeakPoints[i].toInt()] * .70) {
                        break;
                    }
                    if(j == rPeakPoints[i].toInt() + 30) {
                        if(ecgPoints[j]>ecgPoints[j+1])
                            j = rPeakPoints[i].toInt() + 30;
                        else
                            j = rPeakPoints[i].toInt() + 10;
                        break;
                    }
                    j = j + 1;
                }

                while (ecgPoints[j + 1] < ecgPoints[j]) j++
                sPoints.add(j.toDouble())
            }
//        if (isAugmentedLead)
//            sPoints = maxMin(ecgPoints, sPoints, 1)
            return sPoints
        }

        private fun calculateTwaveEndPoints(Tpoints: ArrayList<Double>): ArrayList<Double> {
            val signalToConsider = if (isAugmentedLead) ecgPointsForAugmentedLead!! else ecgPoints

            val tEndPoints = ArrayList<Double>()
            for (i in Tpoints.indices) {
                var j = Tpoints[i].toInt() + 20
                while (signalToConsider[j] / signalToConsider[j + 1] > 1.0001) j++
                tEndPoints.add(j.toDouble())
            }
            return tEndPoints
        }

        private fun calculateQpoints(): ArrayList<Double> {
            val qPoints = ArrayList<Double>()
            for (i in 1 until rPeakPoints.size) {
                var j = rPeakPoints[i].toInt() - 5

                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                if(ecgPoints[rPeakPoints[i].toInt()]>900.0) {
                    for(l in rPeakPoints[i].toInt()  downTo rPeakPoints[i].toInt() -21) {
                        if (ecgPoints[l] < 900.0) {
                            j = l - 5;
                            break
                        }
                    }

                    if(j == rPeakPoints[i].toInt()  - 5) {
                        j = rPeakPoints[i].toInt()  - 25;
                    }

                }
                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


                while (ecgPoints[j] / ecgPoints[j - 1] > 1.0003) //ADDED ON 10-12-2018
                    j--
                qPoints.add(j.toDouble())
            }
//        if (isAugmentedLead)
//            qPoints = maxMin(ecgPoints, qPoints, 1)

            return qPoints
        }

        private fun calculateQRSstopPoints(
            sPoints: ArrayList<Double>,
            prStopPoints: ArrayList<Double>?
        ): ArrayList<Double> {
            val signalToConsider = if (isAugmentedLead) ecgPointsForAugmentedLead!! else ecgPoints

            val qrsStopPoints = ArrayList<Double>()
            for (i in 0 until sPoints.size - 1) {
                var j = sPoints[i + 1].toInt() + 2
                while ((signalToConsider[j + 1] / signalToConsider[j]) > 1.00099 && (signalToConsider[j] <= signalToConsider[prStopPoints!![i].toInt()])) {
                    j++
                }
                /**10-NOV-2021
                 * qrsStopPoints entry increament by 30*/
                qrsStopPoints.add(j.toDouble() + 15)//18 feb 30>15
            }
            return qrsStopPoints
        }

        private fun calculatePPoints(qPoints: ArrayList<Double>): ArrayList<Double> {
            val pPoints = ArrayList<Double>()
            for (i in 0 until rPeakPoints.size - 1) {
                val currentRRInterval = (rPeakPoints[i + 1] - rPeakPoints[i]).toInt()
//            if (leadPosition == 0)
//                println("MIDDLE: ${(rPeakPoints[i] + (0.65 * currentRRInterval).roundToInt())}")
                val pPoint = calculateMaxPoint(
                    (rPeakPoints[i] + (0.65 * currentRRInterval).roundToInt()).toInt(),
                    qPoints[i].toInt()
                ).toDouble()
                pPoints.add(pPoint)
            }
//        if (isAugmentedLead)
//            pPoints = maxMin(ecgPoints, pPoints, 0)

            return pPoints
        }
        private fun calculateTpointsOld(sPoints: ArrayList<Double>): ArrayList<Double> {
//        if(leadPosition==6)
//            println()
            val tPoints = ArrayList<Double>()
            for (i in 0 until rPeakPoints.size - 1) {

                val currentRRInterval = (rPeakPoints[i + 1] - rPeakPoints[i]).toInt()
//            if (leadPosition == 7)
//                println("RR number ")
                /**
                 * date 10-NOV-2021
                 * constant value 0.55 changed to 0.35 */

                val pPoint = calculateMaxPoint(
                    sPoints[i].toInt(),
                    (rPeakPoints[i] + (0.55 * currentRRInterval)).roundToInt()
                ).toDouble()
//            println(pPoint)
                tPoints.add(pPoint)
            }
//        if (isAugmentedLead)
//            tPoints = maxMin(ecgPoints, tPoints, 0)

            return tPoints
        }
//    private fun calculateTpoints(sPoints: ArrayList<Double>): ArrayList<Double> {
//        val tPoints = ArrayList<Double>()
//        for (i in 0 until rPeakPoints.size - 1) {
//
//            val currentRRInterval = (rPeakPoints[i + 1] - rPeakPoints[i]).toInt()
////            if (leadPosition == 7)
////                println("RR number ")
//            /**
//             * date 10-NOV-2021
//             * constant value 0.55 changed to 0.35 */
//
//            val base=ecgPoints[(rPeakPoints[i] + (0.55 * currentRRInterval)).roundToInt()];
//            var j=1
//            var start1=sPoints[i].toInt()+5;
//            var start2=sPoints[i].toInt()-5;
//
//            while(j!=-1)
//            {
//                if( start1+j<rPeakPoints[i+1].toInt())
//                {
//                    var slope1=ecgPoints[start1+j+1]-ecgPoints[start1+j];
//                    var slope2=ecgPoints[start2-j-1]-ecgPoints[start2-j];
//                    if(abs((slope2-slope1)/slope2) >.98 && abs(slope1)<(1/sqrt(3.0)) )
//                    {
//                        start1 += j;
//                        j = -2;
//                    }
//                }
//                else {
//                    j = -2;
//                }
//
//                j += 1;
//            }
//
//            var stop=ceil(rPeakPoints[i] + (0.55 * currentRRInterval)).toInt();
//            var pPoint:Double
//            if(start1<stop) {
//                pPoint = calculateMaxPointTWave(
//                    start1,
//                    stop
//                ).toDouble()
//            }
//            else {
//                pPoint = calculateMaxPoint(
//                    sPoints[i].toInt(),
//                    stop
//                ).toDouble()
//            }
////            println(pPoint)
//            tPoints.add(pPoint)
//        }
////        if (isAugmentedLead)
////            tPoints = maxMin(ecgPoints, tPoints, 0)
//
//        return tPoints
//    }

        private fun calculateTpoints(): ArrayList<Double> {
            var tPoints = ArrayList<Double>()

            var isAnyParameterEmpty=false

            if(rPeakPoints.size==0 || prStopIndices.size==0 || qrsStopIndices.size==0 || ecgPoints.size==0 || sPoints.size==0) {
                isAnyParameterEmpty = true
            }

            if(!isAnyParameterEmpty)
            {
                var basePoints=ArrayList<Double>();
                basePoints.add(0.0);
                for(i in 0 until prStopIndices.size)
                {
                    basePoints.add(prStopIndices[i]);
                }
                basePoints.add(ecgPoints.size.toDouble()-1.0);


                var baseLine=ArrayList<Double>();
                for(i in 0 until basePoints.size-1) {
                    var x1 = basePoints[i]
                    var x2 = basePoints[i + 1]
                    var y1 = ecgPoints[x1.toInt()]
                    if(i==0)
                        y1 = ecgPoints[basePoints[i + 1].toInt()]
                    var y2 = ecgPoints[x2.toInt()]
                    if(i==basePoints.size-1)
                        y2 = ecgPoints[basePoints[i].toInt()]
                    var m1 = 0
                    var m2 = 0
                    var a = y1;
                    var b = m1;
                    var c = 0.0
                    var d = 0.0
                    if (x2 != x1) {
                        c = (3 * (y2 - y1) / Math.pow((x2 - x1).toDouble(), 2.0)) - ((2 * m1 + m2) / (x2 - x1))
                        d = ((m1 + m2) / Math.pow((x2 - x1).toDouble(), 2.0)) - (2 * (y2 - y1) / Math.pow(
                            (x2 - x1).toDouble(),
                            3.0
                        ))
                    }

                    for (j in x1.toInt() until x2.toInt() - 1) {
                        baseLine.add(
                            a + b * (j - x1) + c * Math.pow(
                                (j - x1).toDouble(),
                                2.0
                            ) + d * Math.pow((j - x1).toDouble(), 3.0));
                    }
                }

                baseLine.add(ecgPoints[ecgPoints.size - 1]);

                var first:Double;
                var last:Double;
                for(i in 0 until rPeakPoints.size-1) {
                    if (i == 0) {
//                first=2*sPoints[i]-rPeakPoints[i];
                        var toggle=true
                        first=Math.ceil((sPoints[i]+rPeakPoints[i])/2);
//                if(ecgPoints[sPoints[i].toInt()+10]>ecgPoints[sPoints[i].toInt()]) {
                        for (q in sPoints[i].toInt() until ceil(rPeakPoints[i] + .35 * (rPeakPoints[i + 1] - rPeakPoints[i])).toInt()) {
                            if ((ecgPoints[q] < baseLine[q] && ecgPoints[q + 1] > baseLine[q + 1]) || (ecgPoints[q] == baseLine[q] && ecgPoints[q + 1] > baseLine[q + 1])) {
                                first = q.toDouble();
                                toggle=false
                                break
                            }
                        }
//                }
//                else
//                {
                        if(toggle) {
                            for (q in sPoints[i].toInt() downTo rPeakPoints[i].toInt()) {
                                if ((ecgPoints[q] < baseLine[q] && ecgPoints[q - 1] > baseLine[q - 1]) || (ecgPoints[q] == baseLine[q] && ecgPoints[q - 1] > baseLine[q - 1])) {
                                    first = q.toDouble();
                                    break
                                }
                            }
                        }
//                }
                    }
                    else {
                        first = qrsStopIndices[i - 1];

                        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        var countFraction = 0.0;
                        for(j in  qrsStopIndices[i - 1].toInt() + 1 until qrsStopIndices[i-1].toInt()+11) {
                            if(ecgPoints[j] > ecgPoints[qrsStopIndices[i - 1].toInt()])
                            {
                                countFraction = countFraction + 1.0;
                            }
                        }

                        if(countFraction / 10.0 > .5) {
                            for(q in first.toInt() until ceil(rPeakPoints[i]+.35 * (rPeakPoints[i+1]-rPeakPoints[i])).toInt()+1) {
                                if(ecgPoints[q] < baseLine[q] && ecgPoints[q + 1] > baseLine[q + 1] ||
                                    ecgPoints[q] == baseLine[q] && ecgPoints[q + 1] > baseLine[q + 1]) {
                                    first = q.toDouble();
                                    break
                                }
                            }
                        }
                        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    }

                    last = ceil((prStopIndices[i]+(rPeakPoints[i] + .55 * (rPeakPoints[i + 1] - rPeakPoints[i])))/2);
                    var max1 = 0.0;
//            var carry = sPoints[i];
                    var carry:Double
                    if(sPoints[i]>first){
                        carry = ceil((first+last)/2.0);
                    }
                    else{
                        carry = ceil((sPoints[i]+last)/2.0);
                    }

                    if(first!=0.0&&last!=0.0) {
                        var count = ArrayList<Double>();
                        var sign = -1;
                        for(j in first.toInt() until last.toInt()) {
                            if(j>sPoints[i])
                                count.add(ecgPoints.get(j) - baseLine.get(j));
                        }

                        if(abs(Filters.max(count)) > abs(Filters.min(count)))
                            sign = 1;

                        for(iterate in 0 until 2) {
                            for (j in first.toInt() until last.toInt()) {
                                if (j > sPoints[i]) {
                                    if (abs(ecgPoints.get(j) - baseLine.get(j)) >= max1) {

                                        if (sign * ecgPoints.get(j) > sign * ecgPoints.get(j - 1) && sign * ecgPoints.get(j) > sign * ecgPoints.get(
                                                j + 1
                                            ) ||
                                            sign * ecgPoints.get(j) == sign * ecgPoints.get(j - 1) && sign * ecgPoints.get(j) > sign * ecgPoints.get(
                                                j + 1
                                            ) ||
                                            sign * ecgPoints.get(j) > sign * ecgPoints.get(j - 1) && sign * ecgPoints.get(j) == sign * ecgPoints.get(
                                                j + 1
                                            )
                                        ) {
                                            carry = j.toDouble();
                                            max1 = abs(ecgPoints.get(j) - baseLine.get(j));
                                        }

                                    }
                                }
                            }
                            if (max1 == 0.0) {
                                sign *= -1;
                            }else
                                break;
                        }
                        tPoints.add(carry);
                    }
                }
            }
            else
                tPoints=calculateTpointsOld(sPoints)

//        if(leadPosition==1)
//            println()
            return tPoints
        }


        private fun calculateMaxPoint(start: Int, stop: Int): Int {
            var maxPoint = start
            var max = ecgPoints[maxPoint]
            for (i in start+1..stop) {
                if (ecgPoints[i] >= max) {
                    max = ecgPoints[i]
                    maxPoint = i
                }
            }
            return maxPoint
        }

//    private fun calculateMaxPointTWave(start: Int, stop: Int): Int {
//        var max = 0.0
//        var maxPoint = start
//        val x1=start
//        val x2=stop
//        val y1=ecgPoints[x1]
//        val y2=ecgPoints[x2]
//        val m1=(y2-y1)/(x2-x1)
//        val m2=0.0
//        val a=y1
//        val b=m1
//        val c=(3*(y2-y1)/Math.pow((x2-x1).toDouble(),2.0))-((2*m1+m2)/(x2-x1))
//        val d=((m1+m2)/Math.pow((x2-x1).toDouble(),2.0))-(2*(y2-y1)/Math.pow((x2-x1).toDouble(),3.0))
//        var reference=0.0
//        for (i in start+1..stop) {
//            reference = a + b * (i - x1) + c * Math.pow((i - x1).toDouble(),2.0)+d*Math.pow((i - x1).toDouble(),3.0);
//            if (Math.abs(ecgPoints[i] - reference) >= max) {
//                max = Math.abs(ecgPoints[i] - reference)
//                maxPoint = i
//            }
//
//        }
//        return maxPoint
//    }

        public fun calculatePRstartIndices(pPoints: ArrayList<Double>): ArrayList<Double> {
            val signalToConsider = if (isAugmentedLead) ecgPointsForAugmentedLead!! else ecgPoints
            val prStartIndices = ArrayList<Double>()
            val array = ArrayList<Double>()

            for (i in pPoints.indices) {
                var j = pPoints[i].toInt() - 20

                while (signalToConsider[j] / signalToConsider[j - 1] > 1.0005) {
                    j--
                    if(j==0) {
                        j = pPoints[i].toInt() - 20
                        break
                    }
                }

                prStartIndices.add(j.toDouble())
            }
            return prStartIndices
        }


        public fun calculatePRstopIndices(qPoints: ArrayList<Double>): ArrayList<Double> {
            val signalToConsider = if (isAugmentedLead) ecgPointsForAugmentedLead!! else ecgPoints

            val prStopIndices = ArrayList<Double>()
            for (i in qPoints.indices) {
                var j = qPoints[i].toInt() - 2
                var len = 0
                while ((signalToConsider[j - 1] / signalToConsider[j]) > 1.0010 || len < 3) {
                    j--
                    len++
                    if(j==0) {
                        j = qPoints[i].toInt() - 2
                        break
                    }
                }
                if (j - pPoints[i] < 20) j = (qPoints[i] - 2).toInt()
                prStopIndices.add(j.toDouble())
            }
            return prStopIndices
        }

        fun calculateRPoints(signal: ArrayList<Double>): ArrayList<Double> {
            if(leadPosition==6)
                println()

            ////////////////code to make signal normal///////////////////////////


            var countIntervalAbove900=0
            var ii = 0
            var key = 0
            val limitAmp=950
            while(ii<signal.size-1) {
                if(signal[ii] > limitAmp) {
                    key = 1
                    if(signal[ii + 1] < limitAmp)
                        key = -1
                }
                if(key == -1) {
                    countIntervalAbove900 = countIntervalAbove900 + 1
                    key = 0
                }
                ii++
            }

//        var j=0
            if(countIntervalAbove900<=3 && countIntervalAbove900>=1) {
                key = 0;
                for (i in signal.indices) {
                    if (signal[i] > limitAmp)
                    {
                        if (key == 0 && i > 99)
                        {
                            for (j in i - 100 until i )
                                sampleToIgnore[j] = j;
                        }
                        sampleToIgnore[i] = i;
                        key = 1;
                    }
                    else
                    {
                        if (key == 1 && i < signal.size - 199)
                        {
                            for (j in i until i + 201)
                            {
                                sampleToIgnore[j] = j
                            }
                            key = 0
                        }
                    }
                }
            }


            ////////////////code to make signal normal///////////////////////////

            val derivedSignal = calculateDerivative(signal)

            ////////////////code to make signal normal///////////////////////////
            var derivedSignal1=ArrayList<Double>()
            derivedSignal1.addAll(derivedSignal)

            for(i in sampleToIgnore.indices) {
                if(sampleToIgnore[i]!= 0)
                    derivedSignal1[i] = 0.0
            }

            val maximaPoints = calculateMaximaPoints(derivedSignal1)
            ////////////////code to make signal normal///////////////////////////


//        val maximaPoints = calculateMaximaPoints(derivedSignal)
            val realPeaks = ArrayList<Double>()
            val rPeaks = ArrayList<Double>()
            val firstPeak: Int
            val lastPeak: Int
            if (maximaPoints.size > 0) {
                firstPeak = maximaPoints[0].toInt()
                lastPeak = maximaPoints[maximaPoints.size - 1].toInt()
            } else {
                firstPeak = 0
                lastPeak = ecgPoints.size - 1
            }
            var index = firstPeak + 100
            var maxPoint = -99999.0
            var secondMaxPoint = maxPoint
            for (i in 4 until derivedSignal.size - 1) {
                if(i!=sampleToIgnore[i]) {
//            if (derivedSignal[i + 1] <= derivedSignal[i] && derivedSignal[i - 1] < derivedSignal[i] || derivedSignal[i + 1] < derivedSignal[i] && derivedSignal[i - 1] <= derivedSignal[i]) {
                    if (derivedSignal[i + 1] <= derivedSignal[i] && derivedSignal[i - 1] < derivedSignal[i]) {  //added for New Algorithm : 24 Feb 2022
                        if (derivedSignal[i] > secondMaxPoint) {
                            if (derivedSignal[i] > maxPoint) {
                                secondMaxPoint = maxPoint
                                maxPoint = derivedSignal[i]
                            } else {
                                secondMaxPoint = derivedSignal[i]
                            }
                        }
                    }
                }
            }
            val threshold = 0.45 * secondMaxPoint
            while (index <= lastPeak - 102) {
                if (derivedSignal[index] >= threshold) {
                    while (derivedSignal[index + 1] > derivedSignal[index] && index < derivedSignal.size) index++
                    realPeaks.add(index.toDouble())
                    index += 100
//                index += 130 // making shortest change to resolve mixing of R and T peaks
                }
                index++
            }

            var meanEcgPoints = Filters.mean(ecgPoints)
            var realPeakMagnitudes = ArrayList<Double>()
            for (i in realPeaks) {
                realPeakMagnitudes.add(ecgPoints[i.toInt()])
            }

            var meanRealPeakMagnitudes = Filters.mean(realPeakMagnitudes)
////////////////////////////////rPeak Improvement///////////////////////////////////////////////
//        if(meanEcgPoints < meanRealPeakMagnitudes) {
//            realPeaks.clear()
//            realPeakMagnitudes.clear()
//            var i = firstPeak + 100
//            while(i < lastPeak - 102) {
//                if(abs(derivedSignal[i]) >= threshold) {
//                    while(abs(derivedSignal[i + 1]) > abs(derivedSignal[i])) {
//                        i = i + 1;
//                    }
//                    if(derivedSignal[i] < 0) {
//                        var j = i;
//                        while(j > i - 30) {
//                            if(derivedSignal[j] > 0) {
//                                i = j;
//                                break;
//                            }
//                            j = j - 1;
//                        }
//                    }
//                    realPeaks.add(i.toDouble());
//                    i = i + 100;
//                }
//                i = i + 1;
//            }
//        }
////////////////////////////////rPeak Improvement///////////////////////////////////////////////


            var array=ArrayList<Int>();
            var count=0;
            for(i in 0 until realPeaks.size-1) {
                if(abs(ecgPoints.get(realPeaks.get(i + 1).toInt()) - ecgPoints.get(realPeaks.get(i).toInt())) / (realPeaks.get(i + 1) - realPeaks.get(i)) > sqrt(3.0)) {
                    array.add(i +1-count);
                    count=count+1;
                }
            }

            for(i in 0 until array.size)
                realPeaks.removeAt(array.get(i))

////////////////////////////////rPeak Improvement///////////////////////////////////////////////
//         meanEcgPoints = Filters.mean(ecgPoints)
////         realPeakMagnitudes = ArrayList<Double>()
//        for (i in realPeaks) {
//            realPeakMagnitudes.add(ecgPoints[i.toInt()])
//        }
//        meanRealPeakMagnitudes = Filters.mean(realPeakMagnitudes)
////////////////////////////////rPeak Improvement///////////////////////////////////////////////

            var countPositiveSlope=0
            var countWayTOPositiveSlope=0
            for(i in  realPeaks.indices) {
                if(ecgPoints[realPeaks[i].toInt() - 1] < ecgPoints[realPeaks[i].toInt() + 1]) {
                    countPositiveSlope = countPositiveSlope + 1;

                    var right = 0;
                    var j = realPeaks[i].toInt();
                    while (j != -1 && j < ecgPoints.size) {
                        right = right + 1;
                        if (ecgPoints[j] < ecgPoints[realPeaks[i].toInt()]) {
                            j = -2;
                        }
                        j = j + 1;
                    }

//                var left = 0;
//                j = realPeaks[i].toInt();
//                while (j != -1 && j > 1) {
//                    left = left + 1;
//                    if (ecgPoints[j] > ecgPoints[realPeaks[i].toInt()]) {
//                        j = 0;
//                    }
//                    j = j - 1;
//                }

                    if (right < 30) {
                        countWayTOPositiveSlope = countWayTOPositiveSlope + 1;
                    }

                }
            }

            if(countPositiveSlope==realPeaks.size) {
                if(countWayTOPositiveSlope == countPositiveSlope)
                    countPositiveSlope = 1;
                else
                    countPositiveSlope = 0;
            }

            if (meanEcgPoints < meanRealPeakMagnitudes || countPositiveSlope==1) {
                for (i in realPeaks.indices) {
                    var j = realPeaks[i].toInt()
                    while (derivedSignal[j] > 0 || derivedSignal[j + 1] > 0) j++
                    ////////////////////////////////////////////////////////////////////////////////////////////////////////
                    var carry=j;
//                var l=j;
//                while(l<j+30) {
//                    if(derivedSignal[l] > 0) {
//                        if(ecgPoints[l] > ecgPoints[carry]) {
//                            carry = l;
//                        }
//                    }
//                    l = l + 1;
//                }
                    ////////////////////////////////////////////////////////////////////////////////////////////////////////
                    rPeaks.add(carry.toDouble())

                }
            } else {
                for (i in realPeaks.indices) {
                    var j = realPeaks[i].toInt()
                    while (j != -1  ) {
                        if(j>0 && j<ecgPoints.size) {
                            if (ecgPoints[j] > ecgPoints[j + 1] && ecgPoints[j] > ecgPoints[j - 1] ||
                                ecgPoints[j].equals(ecgPoints[j + 1]) && ecgPoints[j] > ecgPoints[j - 1] ||
                                ecgPoints[j] > ecgPoints[j + 1] && ecgPoints[j].equals(ecgPoints[j - 1])
                            ) {
                                if (ecgPoints[j] > ecgPoints[realPeaks[i].toInt()]) {
                                    rPeaks.add(j.toDouble())

                                    if (ecgPoints[realPeaks[i].toInt()] > ecgPoints[realPeaks[i].toInt() - 1] && ecgPoints[realPeaks[i].toInt() + 1] > meanEcgPoints)
                                        j = -2
                                    else
                                        j = 0
                                }
                            }
                            if (ecgPoints[realPeaks[i].toInt()] > ecgPoints[realPeaks[i].toInt() - 1] && ecgPoints[realPeaks[i].toInt() + 1] > meanEcgPoints)
                                j += 1
                            else
                                j -= 1
                        }
                        else {
                            rPeaks.add(realPeaks[i]);
                            break;
                        }
                    }
                }
            }

            var countPattern=0;
            var i=1;
            while(i<rPeaks.size-1) {
                if ((rPeaks[i] - rPeaks[i - 1]) / (rPeaks[i + 1] - rPeaks[i - 1]) < .4 || (rPeaks[i] - rPeaks[i - 1]) / (rPeaks[i + 1] - rPeaks[i - 1]) > .6) {
                    countPattern = countPattern + 1;
                    i = i + 1;
                }
                i = i + 1;
            }

            if(countPattern.toDouble()/(rPeaks.size-2)>=.3 && countPattern.toDouble()/(rPeaks.size-2)<=.5) {
                var removeRPeak = ArrayList<Int>();
                var countRpeaks = 0;
                for (i in rPeaks.indices) {
                    var j = rPeaks[i].toInt()
                    var left = 0;
                    while (j != -1) {
                        if ((ecgPoints[j] < ecgPoints[j + 1] && ecgPoints[j] < ecgPoints[j - 1]) ||
                            (ecgPoints[j] == ecgPoints[j + 1] && ecgPoints[j] < ecgPoints[j - 1]) ||
                            (ecgPoints[j] < ecgPoints[j + 1] && ecgPoints[j] == ecgPoints[j - 1])
                        ) {
                            left = j;
                            j = 0;
                        }
                        j = j - 1;
                    }

                    j = rPeaks[i].toInt()
                    var right = 0
                    while (j != -1) {
                        if ((ecgPoints[j] < ecgPoints[j + 1] && ecgPoints[j] < ecgPoints[j - 1]) ||
                            (ecgPoints[j] == ecgPoints[j + 1] && ecgPoints[j] < ecgPoints[j - 1]) ||
                            (ecgPoints[j] < ecgPoints[j + 1] && ecgPoints[j] == ecgPoints[j - 1])
                        ) {
                            right = j;
                            j = -2;
                        }
                        j = j + 1;
                    }

                    if (Math.min((rPeaks[i].toInt() - left), right - rPeaks[i].toInt()) * 2  > 60 && ecgPoints[rPeaks[i].toInt()] < 900
                    ) {
                        removeRPeak.add(i - countRpeaks);
                        countRpeaks++;
                    }
                }

                if(removeRPeak.size!=rPeaks.size) {
                    for (i in 0 until removeRPeak.size)
                        rPeaks.removeAt(removeRPeak.get(i))
                }
            }
//        if(leadPosition==1)
//            println()
            return rPeaks

        }

        private fun calculateRRInterval() {
            rrIntervals = ArrayList()
            for (i in 0 until rPeakPoints.size - 1) {
                rrIntervals.add(rPeakPoints[i + 1] - rPeakPoints[i])
            }
        }

        val rrInterval: Int
            get() {
                if (!isSignalProcessed) return 0
                var sum = 0.0
                for (RRInterval in rrIntervals) sum += RRInterval
                if (rrIntervals.size > 0)
                    return (sum / rrIntervals.size * msDifference).roundToInt()
                return 0;
            }


        val firstRPeakPoint: Int
            get() {
                if (!isSignalProcessed) return ecgPoints[0].toInt()
                return if (rPeakPoints.size > 0) rPeakPoints[0].toInt() else 0
            }

        private fun calculateHeartRate() {
            heartRate = if (rrIntervals.size > 0)
                (60000.0 / rrInterval).roundToInt()
            else 0;
        }


        override fun toString(): String {
            return "ECGProcessing{" +
                    ", heartRate=" + heartRate +
                    ", MS_DIFFERENCE=" + msDifference +
                    ", QRS=" + qrsInterval +
                    ", QT=" + qtInterval +
                    ", QTc=" + qTcInterval +
                    ", PR=" + prInterval +
                    ", isSignalProcessed=" + isSignalProcessed +
                    ", STSegmentInMV=" + stSegmentInMv +
                    '}'
        }

        companion object {

            const val RETAKE_MSG = "please retake the test"

            fun isEcgSignalCompatibleForProcessing(ecgPoints: ArrayList<Double>): Boolean {
//            val derivedSignal = calculateDerivative(ecgPoints)
//            val maximaPoints = calculateMaximaPoints(derivedSignal)
//            return maximaPoints.size > 1
                val derivedSignal = calculateDerivative(ecgPoints)
                val maximaPoints = calculateMaximaPoints(derivedSignal)
                var isTheSignalECG=false;
//            var isTheSignalNotInRange=false;
//            for (i in 0 until ecgPoints.size)
//            {
//                if(ecgPoints[i]<995 || ecgPoints[i]>1000)
//                {
//                    isTheSignalNotInRange=true;
//                    break
//                }
//            }

                var high=0.0;
                var low=1000.0;
                for (i in 0 until ecgPoints.size)
                {
                    if(ecgPoints[i]>high)
                    {
                        high=ecgPoints[i];
                    }
                    if(ecgPoints[i]<low)
                    {
                        low=ecgPoints[i];
                    }
                }

                if(maximaPoints.size > 1) {
                    if (high - low > 30)
                        isTheSignalECG = true
                }

                return isTheSignalECG

            }

            private fun calculateMaximaPoints(signal: ArrayList<Double>): ArrayList<Double> {
                val maxima = ArrayList<Double>()
                var maxPoint = Double.MIN_VALUE
                var secondMaxPoint = maxPoint
                for (i in 4 until signal.size - 1) {
                    if (signal[i + 1] <= signal[i] && signal[i - 1] < signal[i] || signal[i + 1] < signal[i] && signal[i - 1] <= signal[i]) {
                        if (signal[i] > secondMaxPoint) {
                            if (signal[i] > maxPoint) {
                                secondMaxPoint = maxPoint
                                maxPoint = signal[i]
                            } else {
                                secondMaxPoint = signal[i]
                            }
                        }
                    }
                }

                if (maxPoint != 0.0) {
                    val threshold = 0.60 * secondMaxPoint
                    var i = 2
                    while (i < signal.size && signal[i] >= threshold) i++
                    while (i < signal.size - 3) {
                        if (signal[i] >= threshold) {
                            while (i < signal.size - 3 && signal[i + 1] > signal[i]) i++ //CHANGE
                            maxima.add(i.toDouble())
                            while (i < signal.size - 3 && signal[i] >= threshold) i++
                        }
                        i++
                    }
                }
                return maxima
            }

            //Function to calculate first derivative of ECG signal.
            private fun calculateDerivative(signal: ArrayList<Double>): ArrayList<Double> {
                val derivative = ArrayList<Double>()
                for (i in 0 until signal.size - 1) {
                    derivative.add(signal[i + 1] - signal[i])
                }
                return derivative
            }

            fun getFirstPeakIndex(ecgPoints: ArrayList<Double>): Int {
                val derivedSignal = calculateDerivative(ecgPoints)
                val maximaPoints = calculateMaximaPoints(derivedSignal)
                val realPeaks = ArrayList<Double>()
                val Rpeaks = ArrayList<Double>()
                val firstPeak: Int
                val lastPeak: Int
                if (maximaPoints.size > 0) {
                    firstPeak = maximaPoints[0].toInt()
                    lastPeak = maximaPoints[maximaPoints.size - 1].toInt()
                } else {
                    firstPeak = 0
                    lastPeak = ecgPoints.size - 1
                }
                var index = firstPeak + 100
                var maxPoint = -99999.0
                var secondMaxPoint = maxPoint
                for (i in 4 until derivedSignal.size - 1) {
                    if (derivedSignal[i + 1] <= derivedSignal[i] && derivedSignal[i - 1] < derivedSignal[i] || derivedSignal[i + 1] < derivedSignal[i] && derivedSignal[i - 1] <= derivedSignal[i]) {
                        if (derivedSignal[i] > secondMaxPoint) {
                            if (derivedSignal[i] > maxPoint) {
                                secondMaxPoint = maxPoint
                                maxPoint = derivedSignal[i]
                            } else {
                                secondMaxPoint = derivedSignal[i]
                            }
                        }
                    }
                }
                val threshold = 0.45 * secondMaxPoint
                while (index <= lastPeak - 102) {
                    if (derivedSignal[index] >= threshold) {
                        while (derivedSignal[index + 1] > derivedSignal[index]) index++
                        realPeaks.add(index.toDouble())
                        index += 100
                    }
                    index++
                }
                for (i in realPeaks.indices) {
                    var j = realPeaks[i].toInt()
                    while (derivedSignal[j] > 0 || derivedSignal[j + 1] > 0) j++
                    Rpeaks.add(j.toDouble())
                }
                return Rpeaks[0].toInt()
            }
        }

//    fun maxMin(signal: ArrayList<Double>, points: ArrayList<Double>, key: Int): ArrayList<Double> {
//        var i = 0
//        var m = 0
//        var j: Int
//        var index1 = 0
//        val takeT = ArrayList<Double>()
//        while (i != -1) {
//            if (i == points[m].toInt()) {
//                m++
//                j = i - 10
//                val array = java.util.ArrayList<Double>()
//                val arrayT = java.util.ArrayList<Int>()
//                while (j != -1) {
//                    array.add(signal[j])
//                    arrayT.add(j)
//                    if (j == i + 10) {
//                        j = -2
//                    }
//                    j++
//                }
//                val pass = DoubleArray(array.size)
//                val passT = IntArray(arrayT.size)
//                for (g in array.indices) {
//                    pass[g] = array[g]
//                    passT[g] = arrayT[g]
//                }
//                takeT.add(sortT(pass, passT, key).toDouble())
//                index1++
//            }
//            i++
//            if (points.size > 0 && i > points[points.size - 1]) {
//                i = -1
//            }
//        }
//        return takeT
//    }

        private fun sortT(Pass: DoubleArray, Pass_T: IntArray, key: Int): Int {
            val c = Pass_T.size
            var temp: Int
            var temp1: Double
            var q: Int
            var r: Int
            q = 0
            while (q < c) {
                r = 0
                while (r < c - 1) {
                    if (Pass[r] < Pass[r + 1]) {
                        temp1 = Pass[r]
                        Pass[r] = Pass[r + 1]
                        Pass[r + 1] = temp1
                        temp = Pass_T[r]
                        Pass_T[r] = Pass_T[r + 1]
                        Pass_T[r + 1] = temp
                    }
                    r++
                }
                q++
            }
            return if (key == 0) {
                Pass_T[0]
            } else {
                Pass_T[Pass_T.size - 1]
            }
        }

//    fun processForAugmentedLead(
//        pPointsLead2: ArrayList<Double>,
//        qPointsLead2: ArrayList<Double>,
//        rPointsLead2: ArrayList<Double>,
//        sPointsLead2: ArrayList<Double>,
//        tPointsLead2: ArrayList<Double>
//    ): EcgCharacteristics {
//
//        this.pPoints = pPointsLead2
//        this.qPoints = qPointsLead2
//        this.rPeakPoints = rPointsLead2
//        this.sPoints = sPointsLead2
//        this.tPoints = tPointsLead2
//
//        return features(
//            rPointsLead2,
//            pPointsLead2,
//            qPointsLead2,
//            sPointsLead2,
//            tPointsLead2
//        )
//
//    }


//    fun features(
//        rPeaks: ArrayList<Double>,
//        pPointIndices: ArrayList<Double>,
//        qPointIndices: ArrayList<Double>,
//        sPointIndices: ArrayList<Double>,
//        tPointIndices: ArrayList<Double>
//    ): EcgCharacteristics {
//        val prStartIndices = calculatePRstartIndices(pPointIndices)
//        val prStopIndices = calculatePRstopIndices(qPointIndices)
//        val qrsStopIndices = calculateQRSstopPoints(sPointIndices, prStopIndices)
//        val tWaveEndIndices = calculateTwaveEndPoints(tPointIndices)
////        val averagePAmplitudeInLead: Double = calculateAverageValues(pPointIndices, prStopIndices,1);
////        val averageQAmplitudeInLead: Double = calculateAverageValues(qPointIndices, prStopIndices,1);
////        val averageSAmplitudeInLead: Double = calculateAverageValues(sPointIndices, prStopIndices,0);
////        val averageTAmplitudeInLead: Double = calculateAverageValues(tPointIndices, prStopIndices,0);
////        val averageRAmplitudeInLead: Double = calculateAverageValues(rPeaks, prStopIndices,0);
//        var totalPR = 0.0
//        var totalQRS = 0.0
//        var totalQT = 0.0
//        var totalRR = 0.0
//        var totalST = 0.0
//        var totalSTE = 0.0
//        val rrIntervals = ArrayList<Double>()
//        var i: Int = 1
//        while (i < rPeaks.size) {
//            totalRR += (rPeaks[i] - rPeaks[i - 1])
//            rrIntervals.add(rPeaks[i] - rPeaks[i - 1])
//            i++
//        }
//        val handle1 = ArrayList<Double>()
//        var g: Int = 0
//        while (g < rrIntervals.size) {
//            handle1.add(rrIntervals[g])
//            g++
//        }
//        i = 0
//        while (i < prStartIndices.size) {
//            totalPR = (totalPR + (prStopIndices[i] - prStartIndices[i]))
//            i++
//        }
//        val qrsIntervalsList = ArrayList<Double>()
//        i = 0
//        while (i < prStartIndices.size - 1) {
//            totalQRS = (totalQRS + ((qrsStopIndices[i] - 15) - prStopIndices[i]))//18 feb 30>15
//            qrsIntervalsList.add((qrsStopIndices[i] - 15) - prStopIndices[i])//18 feb 30>15
//            i++
//        }
//        i = 1
//        while (i < tWaveEndIndices.size) {
//            totalQT = (totalQT + (tWaveEndIndices[i] - prStopIndices[i - 1]))
//            i++
//        }
//        i = 0
//        while (i < tWaveEndIndices.size - 1) {
//            totalST = (totalST + (tWaveEndIndices[i] - qrsStopIndices[i]))
//            i++
//        }
//        i = 0
//        while (i < qrsStopIndices.size) {
//            totalSTE =
//                totalSTE + ecgPoints[qrsStopIndices[i].toInt()] - ecgPoints[prStopIndices[i].toInt()]
//            i++
//        }
//
//
//        val prInterval =
//            if (prStartIndices.size == 0) {
//                0.0
//            } else {
//                totalPR / prStartIndices.size.toDouble() * msDifference
//            }
//        val qrsInterval = if ((prStartIndices.size - 1) == 0) {
//            0.0
//        } else {
//            totalQRS / (prStartIndices.size.toDouble() - 1) * msDifference
//        }
//
//        val qtInterval = if ((prStartIndices.size - 1) == 0) {
//            0.0
//        } else {
//            totalQT / (prStartIndices.size.toDouble() - 1) * msDifference
//        }
//
//        val rrInterval =
//            if ((rPeaks.size - 1) == 0) {
//                0.0
//            } else {
//                totalRR / (rPeaks.size.toDouble() - 1) * msDifference
//            }
////        val stSegment = totalST / tWaveEndIndices.size.toDouble() * msDifference
//
//        val steInUnits = if (qrsStopIndices.size == 0) {
//            0.0
//        } else {
//            totalSTE / qrsStopIndices.size.toDouble()
//        }
//
//        val steInMV = steInUnits * 0.0029296875
//        val qtc: Double
//        val heartRate = 60000 / rrInterval
//        qtc = if (heartRate > 60 && heartRate < 100) {
//            qtInterval / sqrt(rrInterval / 1000.0)
//        } else {
//            qtInterval + 0.154 * (1 - rrInterval / 1000.0) * 1000.0
//        }
//
//        return EcgCharacteristics(
//            prInterval.toInt(),
//            qrsInterval.toInt(),
//            qtInterval.toInt(),
//            qtc,
//            rrInterval.toInt(),
//            heartRate.toInt(),
//            steInMV,
//            qrsIntervalsList,
//            rrIntervals,
//            prStopIndices,
//            prStartIndices,
//            pPointIndices,
//            qPointIndices,
//            sPointIndices,
//            tPointIndices,
//            rPeaks,
//            tWaveEndIndices,
//            averagePAmplitudeInLead,
//            averageQAmplitudeInLead,
//            averageSAmplitudeInLead,
//            averageTAmplitudeInLead,
//            averageRAmplitudeInLead,
//            pWidth,
//            tWidth,
//            qrsDirectionUpward,
//            valueRS,
//            ventricularActivationLOR,
//            ventricularActivationROR,
//            concavity,
//            frequencyOfPatternInQRS,
//            frequencyOfPatternInRR,
//            pAmplitudeArrayInMv,
//            TRRatioSatisfy,
//            TSRatioSatisfy,
//            multiPWave,
//            pPeaksIndex,
//            prIntervalInMilliSec,
//            arrayRateOfMultiPWave,
//            pWavePresent,
//            tentedTWave
//        )
//    }

    }
