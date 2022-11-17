%Experimental Fourier Transform
%By: Sean Perkins
clear all
%==========================================================================
prompt = {'Original Image','Reduced Image','2D Fourier Transform','kx Slice','kx ft inverse','ky Slice','ky ft inverse','misc'};
defaultinput = {'0','1','2','3','0','4','0','0'};
select = inputdlg(prompt,'Select Plots (Include Order)',[1 50],defaultinput);
a = str2double(select(1,1));
b = str2double(select(2,1));
c = str2double(select(3,1));
d = str2double(select(4,1));
e = str2double(select(5,1));
f = str2double(select(6,1));
g = str2double(select(7,1));
h = str2double(select(8,1));
select = [a;b;c;d;e;f;g;h];
array = [a;b;c;d;e;f;g;h];
NumColumns = max(select);
%==========================================================================SelectionMenu
z = 1;
FFTSHIFT = 1;
while z ~= 0
    MENU1 = menu(sprintf('Choose Display Method'),'Plot Numbers Seperately','Plot Numbers Together','Turn Off FFT Shift');
    switch MENU1
        case 1 % Seperately
            MaxNumber = str2double(inputdlg('Maximum Number','Max Number',[1 50],{'3'}));
            NumInstances = str2double(inputdlg('Number of Instances','Instances',[1 50],{'1'}))-1;
            %NumInstances = str2double(NumInstances(1,1))-1;
            skip = 0;
            z = 0;
        case 2 % Together
            MaxNumber = str2double(inputdlg('Maximum Number','Max Number',[1 50],{'3'}));
            NumInstances = 0;
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
            skip = 1;
            z = 0;
        case 3 % FFT Shift
            FFTSHIFT = 0;
    end
end
%==========================================================================ImageProcessing
for i = 0:MaxNumber
    Number = string(i);
    Digit = strcat("/",string(i+30));
    if skip == 0
        figure("Name",strcat("Number ",Number));
    elseif skip == 1
        hold on
    end
    for j = 0:NumInstances
        Instance = string(j);
        NumZeros = 5 - numel(num2str(j));
        MaxZero = ['0','0','0','0','0'];
        ZeroText = MaxZero(1:NumZeros);
        FileName = strcat('digits/',Digit,Digit,"_",ZeroText,string(j),'.png');
        OriginalImage = imread(FileName);
        gray = rgb2gray(OriginalImage);
        %------------------------------------------------------------------
        average = mean(mean(gray));
        for p = 1:size(gray,1)
            for q = 1:size(gray,2)
                if gray(p,q) == 255
                    redv = 255.0 ; 
                else
                    redv = 0.0 ; 
                end
                ReducedImage(p,q) = redv - average ; 
            end
        end
        meanReduced = mean(mean(ReducedImage));%what is this
        %------------------------------------------------------------------
        plotting(OriginalImage,ReducedImage,i,j,skip,MaxNumber,NumColumns,NumInstances,a,b,c,d,e,f,g,h,prompt,FFTSHIFT,array)
    end
end
%==========================================================================FourierTransform
function [f_ft,fiV,fiH,VerticalSlice,HorizontalSlice] = FourierTransform(AlteredImage)
%--------------------------------------------------------------------------
%Fourier Transform
f_ft = fft2(AlteredImage);

%Select a vertical slice of the Fourier Transform
VerticalSlice = zeros(size(f_ft,1),1);
VerticalSlice(:,1) = (f_ft(:,1));

%Apply the inverse fourier transform to the vertical slice
fiV = ifft((VerticalSlice));

%Select a horizontal slice of the Fourier Transform
HorizontalSlice = zeros(1,size(f_ft,1));
HorizontalSlice(1,:) = (f_ft(1,:));

%Apply the inverse fourier transform to the horizontal slice
fiH = ifft((HorizontalSlice));
%--------------------------------------------------------------------------
end
%==========================================================================Plotting
function plotting(OriginalImage,ReducedImage,i,j,skip,MaxNumber,NumColumns,NumInstances,a,b,c,d,e,f,g,h,prompt,FFTSHIFT,array)

[f_ft,fiV,fiH,VerticalSlice,HorizontalSlice] = FourierTransform(ReducedImage);
z = (real(fftshift(VerticalSlice)));

if skip == 0
    RowBaseZero = j*NumColumns;
    NumRows = NumInstances+1;
elseif skip == 1
    RowBaseZero = i*NumColumns;
    NumRows = MaxNumber+1;
end


for i = 1:size(array)
    selection = array(i)
    if 
    
end


%------------------------------------------------------------------
if array(1) ~= 0
    subplot(NumRows,NumColumns,RowBaseZero+a)
    imshow(OriginalImage);
    title(prompt(1));
end
%------------------------------------------------------------------
if array(2) ~= 0
    subplot(NumRows,NumColumns,RowBaseZero+b)
    imshow(ReducedImage);
    title(prompt(2));
end
%------------------------------------------------------------------
if array(3) ~= 0
    subplot(NumRows,NumColumns,RowBaseZero+c)
    if FFTSHIFT == 1
        surfc(real(fftshift(f_ft)))
    elseif FFTSHIFT == 0
        surfc(real(f_ft))
    end
    shading interp % this changes the color shading (i.e. gets rid of the grids lines on the surface)
    title(prompt(3));
end
%------------------------------------------------------------------
if array(4) ~= 0
    subplot(NumRows,NumColumns,RowBaseZero+d)
    %(real(VerticalSlice))
    if FFTSHIFT == 1
        plot(real(fftshift(VerticalSlice)))
    elseif FFTSHIFT == 0
        plot(real(VerticalSlice))
    end
    title(prompt(4));
end
%------------------------------------------------------------------
if array(5) ~= 0
    subplot(NumRows,NumColumns,RowBaseZero+e)
    plot(fiV)
    title(prompt(5));
end
%------------------------------------------------------------------
if array(6) ~= 0
    subplot(NumRows,NumColumns,RowBaseZero+f)
    %plot(real(HorizontalSlice))
    if FFTSHIFT == 1
        plot(real(fftshift(HorizontalSlice)))
    elseif FFTSHIFT == 0
        plot(real(HorizontalSlice))
    end
    title(prompt(6));
end
%------------------------------------------------------------------
if array(7) ~= 0
    subplot(NumRows,NumColumns,RowBaseZero+g)
    plot(fiH)
    title(prompt(7));
end
%------------------------------------------------------------------
if array(8) ~= 0
    subplot(NumRows,NumColumns,RowBaseZero+h)
end
%------------------------------------------------------------------
end
%==========================================================================