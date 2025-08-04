# SARIMA-F-
# SARIMA Implementation in F#

This F# library provides a complete implementation of Seasonal AutoRegressive Integrated Moving Average (SARIMA) models for time series analysis and forecasting. It includes tools for model fitting, parameter optimization, and forecasting with visualization capabilities.

## Features

- **Full SARIMA model implementation**:
  - Non-seasonal and seasonal components
  - Differencing (regular and seasonal)
  - Parameter optimization
- **Model selection**:
  - Automatic order selection with AIC/BIC criteria
- **Forecasting**:
  - Multi-step ahead predictions
  - Visualization of results
- **Synthetic data generation** for testing

## Dependencies

- MathNet.Numerics (v5.0.0+)
- MathNet.Numerics.FSharp
- XPlot.Plotly
- FSharp.Collections.ParallelSeq

## Usage

### 1. Import namespaces
```fsharp
open System
open MathNet.Numerics
open MathNet.Numerics.Optimization
open XPlot.Plotly
```

### 2. Create synthetic time series data
```fsharp
let generateTestSeries n =
    let rng = Random()
    [| for i in 0 .. n - 1 -> 
        let baseVal = sin (2.0 * Math.PI * float i / 20.0)
        let noise = 0.3 * rng.NextDouble()
        baseVal + noise |]

let series = generateTestSeries 200
```

### 3. Fit SARIMA model
```fsharp
let d0, D0, m0 = 1, 1, 5  // Differencing orders and seasonal period
let p, d, q, P, D, Q, model = SARIMA.autoFitSARIMA series d0 D0 m0
```

### 4. Generate forecast
```fsharp
let forecast = SARIMA.forecastSARIMA series model p d q P D Q m0 20
```

### 5. Visualize results
```fsharp
let historicalTrace = Scatter(/* ... */)
let forecastTrace = Scatter(/* ... */)
[historicalTrace; forecastTrace] |> Chart.Plot |> Chart.Show
```

## Key Functions

### SARIMA Module
- `fitSARIMA`: Fit SARIMA model to time series data
- `forecastSARIMA`: Generate forecasts using fitted model
- `autoFitSARIMA`: Automatically select best model orders using AIC
- `combinedDifference`: Apply differencing operations
- `combinedInverseDifference`: Reconstruct original series from differenced data

### Model Parameters
```fsharp
type SARIMAParameters = {
    Mu: float
    Phi: float[]       // Non-seasonal AR coefficients
    Theta: float[]     // Non-seasonal MA coefficients
    SeasonalPhi: float[]   // Seasonal AR coefficients
    SeasonalTheta: float[] // Seasonal MA coefficients
}
```

## Model Selection

The `autoFitSARIMA` function automatically tests combinations of these parameters:
- `p`: Non-seasonal AR order (0-2)
- `d`: Non-seasonal differencing order
- `q`: Non-seasonal MA order (0-2)
- `P`: Seasonal AR order (0-2)
- `D`: Seasonal differencing order
- `Q`: Seasonal MA order (0-2)
- `m`: Seasonal period

It selects the model with the lowest Akaike Information Criterion (AIC).

## Output Example

```
Best model: SARIMA(2,1,2)(2,1,0)[5]
Parameters:
  Mu:    -0.000073
  Phi:   0.775224; 0.100577
  Theta: -1.834099; 0.853107
  Seasonal Phi:   -0.299156; -0.179184
  Seasonal Theta: 

Forecast:
1: 0.035802
2: 0.474528
3: 0.620313
...
```

## Visualization
![SARIMA Forecast](file:///C:/Users/omegam/AppData/Local/Temp/b6695193-ccfa-44a4-964d-0e7623bafeb0.html)

The output plot shows:
- Blue line: Historical data
- Red dashed line: Forecasted values

## Requirements

- .NET 5.0+ or .NET Core 3.1+
- F# 5.0+

## Limitations

1. Model orders (p, q, P, Q) are limited to 0-2 in auto-fit
2. Requires minimum 50 data points
3. Optimization may fail for complex models
4. Seasonal period must be known in advance

## Mathematical Background

SARIMA models are defined as:
```
ARIMA(p,d,q)(P,D,Q)[m]
Where:
  p = non-seasonal AR order
  d = non-seasonal differencing
  q = non-seasonal MA order
  P = seasonal AR order
  D = seasonal differencing
  Q = seasonal MA order
  m = seasonal period
```

The model equation is:
```
(1 - Σφ_i B^i)(1 - ΣΦ_i B^{im}) (1 - B)^d (1 - B^m)^D y_t = 
(1 + Σθ_i B^i)(1 + ΣΘ_i B^{im}) ε_t
```

## References

1. Box, G. E. P., Jenkins, G. M., Reinsel, G. C., & Ljung, G. M. (2015). Time series analysis: forecasting and control. John Wiley & Sons.
2. Hyndman, R. J., & Athanasopoulos, G. (2018). Forecasting: principles and practice. OTexts.
```

This README includes:
1. Key features and dependencies
2. Basic usage examples
3. Core functionality overview
4. Model selection explanation
5. Sample output format
6. Visualization description
7. Requirements and limitations
8. Mathematical background
9. Academic references

The documentation is structured to help new users understand how to use your SARIMA implementation while providing enough technical details for advanced users to understand the underlying methodology.
