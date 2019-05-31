function [results,store_power,fs_array] = lomb_analysis(binplots,time_array,dt,varargin)

warning('This function has been superseded by scargle_analysis.m, which I shall now run for you')
    
[results,store_power,fs_array] = scargle_analysis(binplots,time_array,dt,varargin);