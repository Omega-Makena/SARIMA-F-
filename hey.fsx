#r "nuget: MathNet.Numerics, 5.0.0"
#r "nuget: MathNet.Numerics.FSharp"
#r "nuget: XPlot.Plotly"
#r "nuget: FSharp.Collections.ParallelSeq"

open System
open MathNet.Numerics
open MathNet.Numerics.Optimization
open XPlot.Plotly
open FSharp.Collections.ParallelSeq

module SARIMA =

    type SARIMAParameters = {
        Mu: float
        Phi: float[]
        Theta: float[]
        SeasonalPhi: float[]
        SeasonalTheta: float[]
    }

    let difference (series: float[]) (d: int) = 
        if d = 0 then series
        else
            let rec diff (n: int) (s: float[]) : float[] =
                if n = 0 then s
                else
                    if s.Length < 2 then [||]
                    else
                        s |> Array.pairwise |> Array.map (fun (x1, x2) -> x2 - x1) |> diff (n - 1)
            diff d series

    let seasonalDifference (series: float[]) (D: int) (m: int) =
        if D = 0 then series
        else
            let rec diff (n: int) (s: float[]) : float[] =
                if n = 0 then s
                else
                    if s.Length <= m then [||]
                    else
                        [| for i in m .. s.Length - 1 -> s.[i] - s.[i - m] |] |> diff (n - 1)
            diff D series

    let combinedDifference (series: float[]) (d: int) (D: int) (m: int) : float[] =
        difference (seasonalDifference series D m) d

    let inverseDifference (history: float[]) (diffedForecast: float[]) (d: int) : float[] =
        if d = 0 then diffedForecast
        else
            let result = Array.zeroCreate diffedForecast.Length
            if history.Length < d then
                failwithf "Insufficient history for differencing inversion (need %d, have %d)" d history.Length
            // Start with last known value
            result.[0] <- diffedForecast.[0] + history.[history.Length - d]
            for i in 1 .. diffedForecast.Length - 1 do
                result.[i] <- diffedForecast.[i] + result.[i - 1]
            result

    let inverseSeasonalDifference (history: float[]) (seasonalDiffForecast: float[]) (D: int) (m: int) : float[] =
        if D = 0 then seasonalDiffForecast
        else
            let L = history.Length
            let h = seasonalDiffForecast.Length
            let res = Array.zeroCreate h
            
            let coeffs = 
                [| for k in 1 .. D ->
                    let sign = if k % 2 = 0 then -1.0 else 1.0
                    sign * float (SpecialFunctions.Binomial(D, k)) |]

            for i in 0 .. h - 1 do
                let mutable value = seasonalDiffForecast.[i]
                for k in 1 .. D do
                    let lag = k * m
                    let idx = L + i - lag
                    
                    // Check both history and partial results with bounds
                    if idx >= 0 && idx < L then
                        value <- value + coeffs.[k-1] * history.[idx]
                    elif idx >= L then
                        let resIdx = idx - L
                        if resIdx < i && resIdx >= 0 then
                            value <- value + coeffs.[k-1] * res.[resIdx]
                res.[i] <- value
            res

    let combinedInverseDifference (history: float[]) (forecastDiff: float[]) (d: int) (D: int) (m: int) : float[] =
        // First invert non-seasonal differencing using seasonally differenced history
        let historySeasoned = seasonalDifference history D m
        let invNonSeason = inverseDifference historySeasoned forecastDiff d
        
        // Then invert seasonal differencing using original history
        inverseSeasonalDifference history invNonSeason D m



    let arPart (phi: float[]) (seasonalPhi: float[]) (m: int) (series: float[]) (t: int) : float =
        let arNonSeasonal =
            if not (isNull phi) && phi.Length > 0 then
                phi 
                |> Array.mapi (fun i coeff -> 
                    let idx = t - i - 1
                    if idx >= 0 && idx < series.Length then coeff * series.[idx] else 0.0)
                |> Array.sum
            else 0.0
            
        let arSeasonal =
            if not (isNull seasonalPhi) && seasonalPhi.Length > 0 then
                seasonalPhi 
                |> Array.mapi (fun i coeff -> 
                    let idx = t - m * (i + 1)
                    if idx >= 0 && idx < series.Length then coeff * series.[idx] else 0.0)
                |> Array.sum
            else 0.0
            
        arNonSeasonal + arSeasonal

    let maPart (theta: float[]) (seasonalTheta: float[]) (m: int) (errors: float[]) (t: int) : float =
        let maNonSeasonal =
            if not (isNull theta) && theta.Length > 0 then
                theta 
                |> Array.mapi (fun i coeff -> 
                    let idx = t - i - 1
                    if idx >= 0 && idx < errors.Length then coeff * errors.[idx] else 0.0)
                |> Array.sum
            else 0.0
            
        let maSeasonal =
            if not (isNull seasonalTheta) && seasonalTheta.Length > 0 then
                seasonalTheta 
                |> Array.mapi (fun i coeff -> 
                    let idx = t - m * (i + 1)
                    if idx >= 0 && idx < errors.Length then coeff * errors.[idx] else 0.0)
                |> Array.sum
            else 0.0
            
        maNonSeasonal + maSeasonal

    let residuals (series: float[]) (d: int) (D: int) (m: int)
                  (phi: float[]) (theta: float[]) (seasonalPhi: float[]) (seasonalTheta: float[]) (mu: float) : float[] =
        let diffedSeries = combinedDifference series d D m
        let n = diffedSeries.Length
        let errors = Array.zeroCreate n
        let fitted = Array.zeroCreate n

        let p = if phi <> null then phi.Length else 0
        let q = if theta <> null then theta.Length else 0
        let P = if seasonalPhi <> null then seasonalPhi.Length else 0
        let Q = if seasonalTheta <> null then seasonalTheta.Length else 0
        
        let maxLag = 
            [p; q; P * m; Q * m] 
            |> List.filter (fun x -> x > 0)
            |> fun xs -> if xs.IsEmpty then 0 else List.max xs
        
        if maxLag >= n then 
            failwithf "Time series too short for model order. Required: %d, Available: %d" maxLag n

        // Initialize errors to zero for first maxLag points
        for i in 0 .. maxLag - 1 do
            errors.[i] <- 0.0
            fitted.[i] <- diffedSeries.[i]

        for t in maxLag .. n - 1 do
            let arVal = arPart phi seasonalPhi m diffedSeries t
            let maVal = maPart theta seasonalTheta m errors t
            fitted.[t] <- mu + arVal + maVal
            errors.[t] <- diffedSeries.[t] - fitted.[t]

        errors

    let objective (series: float[]) (d: int) (D: int) (m: int) (p: int) (q: int) (P: int) (Q: int) (parameters: float[]) : float =
        let mu = parameters.[0]
        let phi = if p > 0 then parameters.[1..1+p-1] else [||]
        let theta = if q > 0 then parameters.[1+p..1+p+q-1] else [||]
        let seasonalPhi = if P > 0 then parameters.[1+p+q..1+p+q+P-1] else [||]
        let seasonalTheta = if Q > 0 then parameters.[1+p+q+P..1+p+q+P+Q-1] else [||]
        
        try
            let res = residuals series d D m phi theta seasonalPhi seasonalTheta mu
            res |> Array.sumBy (fun e -> e * e)
        with ex -> 
            printfn "Objective failed: %s" ex.Message
            Double.MaxValue

    let fitSARIMA (series: float[]) (p: int) (d: int) (q: int) (P: int) (D: int) (Q: int) (m: int) =
        if series.Length < 50 then 
            failwith "Insufficient data length for SARIMA (min 50 points)"

        let nParams = 1 + p + q + P + Q
        let initialGuess = Array.create nParams 0.1
        initialGuess.[0] <- Array.average series

        let objectiveFunc (v: MathNet.Numerics.LinearAlgebra.Vector<float>) : float =
            objective series d D m p q P Q (v.ToArray())

        let optimizer = new NelderMeadSimplex(1e-6, 10000)
        let initialPoint = MathNet.Numerics.LinearAlgebra.Vector.Build.DenseOfArray initialGuess
        
        try
            let result = optimizer.FindMinimum(ObjectiveFunction.Value(objectiveFunc), initialPoint)
            let paramVector = result.MinimizingPoint.ToArray()
            {
                Mu = paramVector.[0]
                Phi = if p > 0 then paramVector.[1..p] else [||]
                Theta = if q > 0 then paramVector.[p+1..p+q] else [||]
                SeasonalPhi = if P > 0 then paramVector.[p+q+1..p+q+P] else [||]
                SeasonalTheta = if Q > 0 then paramVector.[p+q+P+1..p+q+P+Q] else [||]
            }, result.FunctionInfoAtMinimum.Value
        with ex ->
            printfn "Optimization failed: %s" ex.Message
            reraise()

    let forecastSARIMA (series: float[]) (parameters: SARIMAParameters) (p: int) (d: int) (q: int) (P: int) (D: int) (Q: int) (m: int) (h: int) : float[] =
        let diffed = combinedDifference series d D m
        let n = diffed.Length
        
        // Calculate maximum required lags
        let maxLag = 
            [p; q; P * m; Q * m] 
            |> List.filter (fun x -> x > 0)
            |> fun xs -> if xs.IsEmpty then 0 else List.max xs
        
        if n <= maxLag then
            failwithf "Differenced series too short for forecasting. Required: %d, Available: %d" (maxLag+1) n
        
        let forecastErrors = Array.zeroCreate (n + h)
        let forecasted = Array.append diffed (Array.zeroCreate h)
        
        // Calculate known errors
        try
            let knownErrors = residuals series d D m parameters.Phi parameters.Theta parameters.SeasonalPhi parameters.SeasonalTheta parameters.Mu
            Array.blit knownErrors 0 forecastErrors 0 (min knownErrors.Length n)
        with ex -> 
            printfn "Residual calculation failed: %s" ex.Message
        
        // Forecast future values
        for t in n .. n + h - 1 do
            let arVal = arPart parameters.Phi parameters.SeasonalPhi m forecasted t
            let maVal = maPart parameters.Theta parameters.SeasonalTheta m forecastErrors t
            forecasted.[t] <- parameters.Mu + arVal + maVal
        
        // Reconstruct forecast
        combinedInverseDifference series forecasted.[n..] d D m

    let aic (rss: float) (k: int) (n: int) : float = 
        let nf = float n
        float (2 * k) + nf * Math.Log(rss / nf)

    let bic (rss: float) (k: int) (n: int) : float = 
        let nf = float n
        float k * Math.Log(nf) + nf * Math.Log(rss / nf)

    let autoFitSARIMA (series: float[]) (d: int) (D: int) (m: int) =
        let orders = [0..2] 
        let mutable bestModel = None
        let mutable bestAIC = Double.PositiveInfinity
        
        printfn "Starting SARIMA model search..."
        printfn "d=%d, D=%d, m=%d" d D m
        
        orders
        |> List.map (fun p -> [p], orders)
        |> List.collect (fun (ps, qs) -> qs |> List.map (fun q -> ps@[q], orders))
        |> List.collect (fun (pqs, Ps) -> Ps |> List.map (fun P -> pqs@[P], orders))
        |> List.collect (fun (pqPs, Qs) -> Qs |> List.map (fun Q -> pqPs@[Q]))
        |> List.map (function 
            | [p; q; P; Q] -> (p, q, P, Q)
            | _ -> failwith "Invalid order combination")
        |> List.filter (fun (p, q, P, Q) -> 
            let totalParams = 1 + p + q + P + Q
            let minDataNeeded = totalParams + d + D * m + 10
            series.Length > minDataNeeded)
        |> List.iter (fun (p, q, P, Q) ->
            try
                printfn "Trying SARIMA(%d,%d,%d)(%d,%d,%d)[%d]" p d q P D Q m
                let model, rss = fitSARIMA series p d q P D Q m
                let k = 1 + p + q + P + Q
                let n = series.Length - d - D * m
                let aicVal = aic rss k n
                let bicVal = bic rss k n
                printfn "SARIMA(%d,%d,%d)(%d,%d,%d)[%d] -> RSS: %.4f, AIC: %.4f, BIC: %.4f" p d q P D Q m rss aicVal bicVal
                
                if aicVal < bestAIC then
                    bestAIC <- aicVal
                    bestModel <- Some (p, d, q, P, D, Q, model)
            with ex ->
                printfn "Model SARIMA(%d,%d,%d)(%d,%d,%d)[%d] failed: %s" p d q P D Q m ex.Message)

        match bestModel with
        | Some model -> model
        | None -> failwith "No SARIMA model converged"


// Test with synthetic data
let generateTestSeries n =
    let rng = Random()
    [| for i in 0 .. n - 1 -> 
        let baseVal = sin (2.0 * Math.PI * float i / 20.0)
        let noise = 0.3 * rng.NextDouble()
        baseVal + noise |]

let n = 200
let series = generateTestSeries n
let d0, D0, m0 = 1, 1, 5  // or m0 = your real seasonal period

printfn "Starting SARIMA fitting..."
let p, d, q, P, D, Q, model = SARIMA.autoFitSARIMA series d0 D0 m0

printfn "\nBest model: SARIMA(%d,%d,%d)(%d,%d,%d)[%d]" p d q P D Q m0
printfn "Parameters:"
printfn "  Mu:    %.6f" model.Mu
printfn "  Phi:   %s" (model.Phi |> Array.map (sprintf "%.6f") |> String.concat "; ")
printfn "  Theta: %s" (model.Theta |> Array.map (sprintf "%.6f") |> String.concat "; ")
printfn "  Seasonal Phi:   %s" (model.SeasonalPhi |> Array.map (sprintf "%.6f") |> String.concat "; ")
printfn "  Seasonal Theta: %s" (model.SeasonalTheta |> Array.map (sprintf "%.6f") |> String.concat "; ")

let forecast = SARIMA.forecastSARIMA series model p d q P D Q m0 20

printfn "\nForecast:"
forecast |> Array.iteri (fun i f -> printfn "%d: %.6f" (i+1) f)

// Plot results
let historicalTrace =
    Scatter(
        x = [| for i in 0 .. series.Length - 1 -> float i |],
        y = series,
        mode = "lines",
        name = "Historical",
        line = Line(color = "blue", width = 2.0)
    )

let forecastTrace =
    Scatter(
        x = [| for i in series.Length .. series.Length + forecast.Length - 1 -> float i |],
        y = forecast,
        mode = "lines+markers",
        name = "Forecast",
        line = Line(color = "red", width = 2.0, dash = "dash")
    )

let chart =
    [historicalTrace; forecastTrace]
    |> Chart.Plot
    |> Chart.WithTitle (sprintf "SARIMA(%d,%d,%d)(%d,%d,%d)[%d] Forecast" p d q P D Q m0)
    |> Chart.WithXTitle "Time"
    |> Chart.WithYTitle "Value"

chart.Show()

