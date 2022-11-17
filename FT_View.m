%Fourier Transform - Batch Conversion
%By: Sean Perkins
clear all
%==========================================================================VariableInput
MaxNumber = 9;
NumInstances = 100;
DesiredMaxAmplitude = 100;%%%%%%%%%%%%%%%%SCALE TO DESIRED AMPLITUDE
%==========================================================================DesignateOutput
InputFolder = "digit";
OutputFolder = "Output/";
DataBatch = "1";
%==========================================================================CreateDirectories
DataPath = strcat(OutputFolder,"Data_Collection_",DataBatch,"/");mkdir(DataPath);
SubDataPath = strcat(DataPath);mkdir(SubDataPath);
ImagePath = strcat(SubDataPath,"_Images/");mkdir(ImagePath);
CSVPath = strcat(SubDataPath,"CSV/");mkdir(CSVPath);
%==========================================================================Main
for i = 0:MaxNumber
    Number = string(i);
    Digit = strcat("/",string(i+30));
    VkxData = zeros(128,NumInstances);
    HkyData = zeros(128,NumInstances);
    %----------------------------------------------------------------------
%     SubDataPath = strcat(DataPath,Number,"/");
%     mkdir(SubDataPath);
%     InstImagePath = strcat(SubDataPath,Number,"_Images/");
%     mkdir(InstImagePath);
    InstImagePath = strcat(ImagePath,Number);
    mkdir(InstImagePath);
    %----------------------------------------------------------------------
    for j = 0:NumInstances
        Instance = string(j);
        NumZeros = 5 - numel(num2str(j));
        MaxZero = ['0','0','0','0','0'];
        ZeroText = MaxZero(1:NumZeros);
        FileName = strcat(InputFolder,Digit,Digit,"_",ZeroText,string(j),'.png');
        clear NumZeros MaxZero ZeroText
        PicName = strcat("/",Number,"-",Instance);
        %------------------------------------------------------------------
        OriginalImage = imread(FileName);
        gray = rgb2gray(OriginalImage);
        average = mean(mean(gray));
        for p = 1:size(gray,1)
            for q = 1:size(gray,2)
                if gray(p,q) == 255
                    redv = 255.0 ;
                else
                    redv = 0.0 ; 
                end
                Reduced(p,q) = redv - average ;
            end
        end
        %meanReduced = mean(mean(Reduced));
        [Vkx,Hky] = FourierTransform(Reduced);
        VkxData(:,j+1) = Vkx;
        HkyData(:,j+1) = Hky;
        imwrite(gray,strcat(InstImagePath,PicName,"_Image.png"));
    end
    %----------------------------------------------------------------------
    [kx_PosData] = DataScaling(VkxData,DesiredMaxAmplitude);
    [ky_PosData] = DataScaling(HkyData,DesiredMaxAmplitude);
    %Excel File = .xlsx, Matlab File = .mat"
    writematrix(kx_PosData, strcat(CSVPath,Number,"_kx.csv"));
    writematrix(ky_PosData, strcat(CSVPath,Number,"_ky.csv"));
end


S = dir(fullfile(CSVPath,'*.CSV'));
for k = 1:numel(S)
    [~,~,mat] = xlsread(fullfile(CSVPath,S(k).name));
    [~,fnm] = fileparts(S(k).name);
    xlswrite(strcat(SubDataPath,'kxkyData.xlsx'), mat, fnm);
end

%==========================================================================FourierTransform
function [Vkx,Hky] = FourierTransform(Reduced)
%--------------------------------------------------------------------------
%Fourier Transform
f_ft = fft2(Reduced);
%Select a vertical slice of the Fourier Transform
VerticalSlice = zeros(size(f_ft,1),1);
VerticalSlice(:,1) = (f_ft(:,1));
%Select a horizontal slice of the Fourier Transform
HorizontalSlice = zeros(1,size(f_ft,1));
HorizontalSlice(1,:) = (f_ft(1,:));
%Adjust/Rename Slice Data
Vkx = real(fftshift(VerticalSlice));
Hky = real(fftshift(HorizontalSlice));
end
%==========================================================================
function [PosData] = DataScaling(DataIn,DesiredMaxAmplitude)
MaxAmp = max(max(DataIn));
ScaleMult = (DesiredMaxAmplitude)/MaxAmp;
PosData = zeros(size(DataIn));
for j = 1:size(DataIn,2)
    for i = 1:size(DataIn,1)
        Pos = round(DataIn(i,j)*ScaleMult);
        PosData(i,j) = Pos;
    end
end
end
%==========================================================================


