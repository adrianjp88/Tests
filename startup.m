clear all;
close all;

win32 = 0;
win64 = 0;

%% paths to 32 bit mex files and dlls to set by user
%-------------------------------------------------
win32 = 1;
GpuFit_path_32 = 'D:\Sources\MexFiles\32';
CpuFit_path_32 = 'D:\Sources\MexFiles\32';
GpuFit_profile_path_32 = 'D:\Sources\MexFiles\32';
CpuFit_profile_path_32 = 'D:\Sources\MexFiles\32';
cminapck_path_32 = 'D:\Sources\MexFiles\32';
GPULMFit_path_32 = 'D:\Sources\MexFiles\32';
%-------------------------------------------------

%% paths to 64 bit mex files and dlls to set by user
%-------------------------------------------------
% win64 = 1;
% GpuFit_path_64 = 'D:\Sources\MexFiles\64';
% CpuFit_path_64 = 'D:\Sources\MexFiles\64';
% GpuFit_profile_path_64 = 'D:\Sources\MexFiles\64';
% CpuFit_profile_path_64 = 'D:\Sources\MexFiles\64';
% cminapck_path_64 = 'D:\Sources\MexFiles\64';
% GPULMFit_path_64 = ''; %non-existent
%-------------------------------------------------


%% add paths
if win32
    addpath (GpuFit_path_32);
    addpath (CpuFit_path_32)
    addpath (GpuFit_profile_path_32)
    addpath (CpuFit_profile_path_32)
    addpath (cminapck_path_32)
    addpath (GPULMFit_path_32)
end
if win64
    addpath (GpuFit_path_64);
    addpath (CpuFit_path_64)
    addpath (GpuFit_profile_path_64)
    addpath (CpuFit_profile_path_64)
    addpath (cminapck_path_64)
    addpath (GPULMFit_path_64) %non-existent
end

addpath ( genpath ( pwd ) )