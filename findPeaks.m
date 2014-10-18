function [peak_x,peak_y] = findPeaks(dataIn,varargin)
% function [peak_x,peak_y] = findPeaks(dataIn,bounds,smoothing,verbose)
% A function to find the peaks in a data set.
% Inputs:
%     dataIn = (nx2 Matrix) Column 1 is the independent variable (x or t).   
%                           Column 2 is the dependent variable (y or f).
%     bounds = Optional - (2x1 vector) Column 1 is the lower bound of the 
%                                      independent variable. Column 2 is 
%                                      the upper bound.
%     smoothing = Optional - (scalar) Number of points to smooth the signal
%                                     if ripple or noise is on signal (you
%                                     may want this if you have too many
%                                     peaks reported back). Anything above
%                                     11 is not recommended for most 
%                                     signals. An odd number will always be
%                                     used - input will be corrected.
%     verbose = Optional - (scalar) 1 or 0 if the peaks are to be plotted
    verbose=0;
    x = dataIn(:,1);
    y = dataIn(:,2);
    bounds = [];
    smoothing_size = [];
    peak_x = [];
    peak_y = [];
    if length(varargin)==3
        bounds = varargin{1};
        smoothing_size = varargin{2};
        verbose = varargin{3};
    elseif length(varargin)==2
        bounds = varargin{1};
        smoothing_size = varargin{2};
    elseif length(varargin)==1
        bounds = varargin{1};
    end
    if ~isempty(bounds)
        lower_idx = find(x>=bounds(1),1,'first');
        upper_idx = find(x<=bounds(2),1,'last');
    else
        lower_idx = 1;
        upper_idx = length(x);
    end
    X = x(lower_idx:upper_idx);
    Y = y(lower_idx:upper_idx);
    
    if ~isempty(smoothing_size)
        smoothing_size = floor(smoothing_size/2);
        y_smooth = smoothing(Y,smoothing_size);
        dydx = derivative(X,y_smooth);
    else
        smoothing_size = 2;
        dydx = derivative(X,Y);
    end
    
    for j=smoothing_size+2:length(dydx)-smoothing_size
        if(dydx(j-1)>0 && dydx(j)<0)
            % if the previous point is decreasing, but this one is
            % increasing....            
            % find the max in the 5-points
            pMax = max(Y(j-smoothing_size:j+smoothing_size));
            A = Y>=pMax;
            B = X>=X(j-smoothing_size);
            C = X<=X(j+smoothing_size);
            % find the index of the max
            pLoc = find(A.*B.*C,1,'first');
            % set the peak values
            peak_x = [peak_x;X(pLoc)]; %#ok<*AGROW>
            peak_y = [peak_y;Y(pLoc)];
        end
    end
    A = [peak_x,peak_y];
    A = unique(A,'rows');
    peak_x = A(:,1);
    peak_y = A(:,2);
    if(verbose==1)
%         figure(100);
        plot(x,y,'b');
        hold on;
        plot(peak_x,peak_y,'r.');
        drawnow;
        hold off;
    end
    
    
    function output = derivative(x,y)
        % 4th order approximation of derivative, central differences
        % 5-point stencil
        n = length(x);
        output = zeros(n,1);
        for i=3:n-2
            numerator = y(i-2)-8*y(i-1)+8*y(i+1)-y(i+2);
            denominator = abs(mean(x(i-2:i+2)))*12;
            output(i) = numerator/denominator;
        end
    end

    function output = smoothing(y,smoove_factor)
        % reduces error from signal noise
        n = length(y);
        output = y;
        for i=smoove_factor+1:n-smoove_factor
            output(i) = mean(y(i-smoove_factor:i+smoove_factor));
        end
    end
            
end