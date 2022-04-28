# DynamicTransportSurrogates
This repository is the code for my final project for CSE 5835. In this work, I implemented two approaches for constructing a surrogate model for a mechanistic model of fluid and solute transport in tumors. The overall goal was to develop the most accurate mapping between biological parameters of interest and an accumulation profile over a fixed time horizon. The first surrogate modeling approach was to learn a model for multidimensional regression between the parameters and the full accumulation profile. The second approach was to learn a mapping between the parameters, the current time step, the current concentration value and the concentration at the next time step. This approach explored the use of multilayer perceptrons (MLP), vanilla recurrent neural networks (RNN), and long short-term memory (LSTM) models.

This project was done in Julia v1.6.1 using Flux.jl v0.12.1,Surrogates v2.2.0, and DifferentialEquations v6.18.0.

There are two files for model assessment in which grid search was used to assess candidate models and determine the appropriate model hyperparameters. The multidimensional regression, full integration model file is BaselineModelGridSearch.ipynb and the dynamic surrogate model file is DynamicModelGridSearch.ipynb. The final model training for each model was done in BaselineFinalModel.ipynb and DynamicFinalModels.ipynb. The full dataset is located in MLFinalData_v3. The final model files are saved as m_baseline.bson, m_MLP_f.bson, m_RNN_f.bson, and m_LSTM_f.bson. These models can be loaded into DataAnalysis.ipynb for plotting figures and calculating performance metrics.

For interested readers follow these instructions to build the final ML models in this project and analyze the results:
1. Download this repository to your local machine
2. Download Julia v1.6.1
3. Launch the Julia REPL, hit the ],and type the following command:
add Flux@v0.12.1,Statistics,Plots,MLDataUtils,DelimitedFiles,StatsPlots,IJulia
4. Hit backspace and run the following commands in the repl
using IJulia; notebook()
5. Navigate to the downloaded repository
6. Launch BaselineFinalModel.ipynb and run all cells ( make sure to uncomment the @save command in the last cell at the very end)
7. Launch DynamicFinalModels.ipynb and run all cells ( make sure to uncomment the @save commands in the last cell at the very end)
8. Launch DataAnalysis.ipynb and run all cells
