data {
    // Dimensional parameters
    int<lower=1> J1; // distinct number of replicates for standard candle.
    int<lower=1> J2; // distinct number of flow rates
    int<lower=1> N1; // Total number of standard candle measurements.
    int<lower=1> N2; // Total number of shock rate measurements.
    int<lower=1, upper=J1> repl[N1]; //  Vector of trial identifier for standard candle
    int<lower=1, upper=J2> shock[N2]; // Vector of trial identifier for shock rate experiments.
}