module RNASeqTests

open BioFSharp.Stats
open Expecto


//let testSeq = seq { for i in 1000. .. 1019. -> ("stringtest"+ i.ToString(),(100.,i))}
let testSeq = seq { for i in 1. .. 2. -> ("stringtest"+ i.ToString(),(i,i))}

let resultRPKM = seq { ("stringtest1", 333333333.3333333); ("stringtest2", 333333333.3333333)}
let resultTPM = seq {("stringtest1", 500000.); ("stringtest2", 500000.)}

[<Tests>]
let RNASeqTests =
    
    testList "RNASeqTests" [
        testCase "RPKM" (fun _ ->
            Expect.equal 
                (RNASeq.rpkms testSeq
                |> Seq.map (fun (x,y) -> x, System.Math.Round (y,4))|> Array.ofSeq)
                (resultRPKM
                |> Seq.map (fun (x,y) -> x, System.Math.Round (y,4)) |> Array.ofSeq)
                "RPKM did not return correct Sequence"
        )
        testCase "TPM" (fun _ ->
            Expect.equal 
                (RNASeq.tpms testSeq 
                |> Seq.map (fun (x,y) -> x, System.Math.Round (y,4))|> Array.ofSeq)
                (resultTPM
                |> Seq.map (fun (x,y) -> x, System.Math.Round (y,4))|> Array.ofSeq)
                "TPM did not return correct Sequence"
        )        
    ]