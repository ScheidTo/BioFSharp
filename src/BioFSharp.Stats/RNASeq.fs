namespace BioFSharp.Stats

open System
open System.Collections.Generic


module RNASeq = 
    
    type RNASeqInput = { 
    GeneID : string
    GeneLength : float
    GeneCount : float
    } with static member Create id gl gc = {GeneID=id;GeneLength=gl;GeneCount=gc}

    type NormalizationMethod = 
        | RPKM
        | TPM

    type NormalizedCounts = { 
    GeneID : string
    NormalizedCount : float
    NormalizationMethod: NormalizationMethod
    } with static member Create id nc nm = {GeneID=id;NormalizedCount=nc;NormalizationMethod=nm}

    let private calcRPM sumOfAllReadsPerMil counts =
        (counts |> float) / sumOfAllReadsPerMil

    let private calcRPKM geneLength rpm = 
        (float rpm) / ((float geneLength) / 1000.) 

    let private rpkmsOf (geneIDs:seq<string>) (length:seq<float>) (counts:seq<float>) =
        let sumOfAllReads = 
                counts
                |> Seq.sum
        let sumOfAllReadsPerMil =
            sumOfAllReads / 1000000. 
        let rpms =
            Seq.map (fun counts -> calcRPM sumOfAllReadsPerMil counts) counts
        let rpkms =
            let rpkm =
                Seq.zip length rpms
                |> Seq.map (fun (length, rpm) -> calcRPKM length rpm)
            rpkm
        let rpkmResult =
            Seq.map2 (fun ids counts -> {GeneID=ids; NormalizedCount=counts; NormalizationMethod=RPKM}) geneIDs rpkms
        rpkmResult

    let rpkms (idLengthAndCounts:seq<RNASeqInput>) =
        rpkmsOf (idLengthAndCounts |> Seq.map (fun x -> x.GeneID))  (idLengthAndCounts |> Seq.map (fun x -> x.GeneLength))  (idLengthAndCounts |> Seq.map (fun x -> x.GeneCount)) 

    let private tpmsOf (idLengthAndCounts:seq<RNASeqInput>) =
        let rpk = 
            idLengthAndCounts
            |> Seq.map (fun idLengthAndCounts -> idLengthAndCounts.GeneCount/idLengthAndCounts.GeneLength/1000.)
        let sumOfAllReads =
            rpk
            |> Seq.sum
        let sumOfAllReadsPerMil =
            sumOfAllReads / 1000000.
        let tpms =
            rpk
            |> Seq.map (fun rpks -> rpks/sumOfAllReadsPerMil)
        let geneID =
            idLengthAndCounts
            |> Seq.map (fun idLengthAndCounts -> idLengthAndCounts.GeneID)
        let tpmResult =
            Seq.map2 (fun ids counts -> {GeneID=ids; NormalizedCount=counts; NormalizationMethod=TPM}) geneID tpms
        tpmResult
    
    let tpms (idLengthAndCounts:seq<RNASeqInput>) =
        tpmsOf idLengthAndCounts 

