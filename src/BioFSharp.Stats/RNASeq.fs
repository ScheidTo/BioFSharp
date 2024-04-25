namespace BioFSharp.Stats

open System
open System.Collections.Generic


module RNASeq = 
    type RNASeqInput = { // type construction
    GeneID : string
    GeneLength : float
    GeneCount : float
    } with static member Create id gl gc = {GeneID=id;GeneLength=gl;GeneCount=gc}
    //module Rpkm =
    let calcRPM sumOfAllReadsPerMil counts =
        (counts |> float) / sumOfAllReadsPerMil

    let calcRPKM geneLength rpm = 
        (float rpm) / ((float geneLength) / 1000.) 

    let rpkms (geneLengthAndCounts:seq<RNASeqInput>) =
        let sumOfAllReads = 
                geneLengthAndCounts
                |> Seq.map (fun geneLengthAndCounts -> geneLengthAndCounts.GeneCount)
                |> Seq.sum
        let sumOfAllReadsPerMil =
            sumOfAllReads / 1000000. 
        let rpms =
            geneLengthAndCounts
            |> Seq.map (fun geneLengthAndCounts -> calcRPM sumOfAllReadsPerMil geneLengthAndCounts.GeneCount)
        let rpkms =
            let geneLengthsRPM = 
                geneLengthAndCounts
                |> Seq.map (fun getLength -> getLength.GeneLength)
            let rpkm =
                Seq.zip geneLengthsRPM rpms
                |> Seq.map (fun (length, rpm) -> calcRPKM length rpm)
            rpkm
        let geneID =
            geneLengthAndCounts
            |> Seq.map (fun geneLengthAndCounts -> geneLengthAndCounts.GeneID)
        let rpkmResult =
            Seq.zip geneID rpkms
        rpkmResult
    //module Tpm =
    let tpms (geneLengthAndCounts:seq<RNASeqInput>) =
        let rpk = 
            geneLengthAndCounts
            |> Seq.map (fun geneLengthAndCounts -> geneLengthAndCounts.GeneCount/geneLengthAndCounts.GeneLength/1000.)
        let sumOfAllReads =
            rpk
            |> Seq.sum
        let sumOfAllReadsPerMil =
            sumOfAllReads / 1000000.
        let tpms =
            rpk
            |> Seq.map (fun rpks -> rpks/sumOfAllReadsPerMil)
        let geneID =
            geneLengthAndCounts
            |> Seq.map (fun geneLengthAndCounts -> geneLengthAndCounts.GeneID)
        let tpmResult =
            Seq.zip geneID tpms
        tpmResult