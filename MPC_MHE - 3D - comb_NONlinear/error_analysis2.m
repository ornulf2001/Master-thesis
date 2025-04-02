clc,clear
addpath("data_20250327_053323")
load("sims.mat")
load("RMSEtotal.mat")
load("RMSEmean.mat")
load("ests.mat")

errorData = sims - ests;

nStates = size