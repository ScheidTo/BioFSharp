namespace BioFSharp.Stats

open System
open System.Collections.Generic


module RNASeq = 
    //module Rpkm =
    let calcRPM sumOfAllReadsPerMil counts =
        (counts |> float) / sumOfAllReadsPerMil

    let calcRPKM geneLength rpm = 
        (float rpm) / ((float geneLength) / 1000.) 

    let rpkms (geneLengthAndCounts:seq<string*(float*float)>) =
        let sumOfAllReads = 
            geneLengthAndCounts
            |> Seq.map (fun (x,(y,z)) -> z)
            |> Seq.sum
        let sumOfAllReadsPerMil =
            sumOfAllReads / 1000000. 
        let rpms = 
            geneLengthAndCounts
            |> Seq.map (fun (head,(l,c)) -> head, (l, calcRPM sumOfAllReadsPerMil c))
        let rpkms = 
            rpms //geneLengthAndCounts // should be rpms i think
            |> Seq.map (fun (head,(l,c)) -> head, calcRPKM l c)
        rpkms

    //module Tpm =
    let tpms (geneLengthAndCounts:seq<string*(float*float)>) =
        let rpk = 
            geneLengthAndCounts
            |> Seq.map (fun (head,(l,c)) -> head, c/l/1000.)
        let sumOfAllReads =
            rpk
            |> Seq.map snd
            |> Seq.sum
        let sumOfAllReadsPerMil =
            sumOfAllReads / 1000000. 
        let tpms =
            rpk
            |> Seq.map (fun (header,rpks) -> header, rpks/sumOfAllReadsPerMil)
        tpms