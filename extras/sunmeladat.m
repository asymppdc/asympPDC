function u = sunmeladat(selection)
%SUNMELADAT   Return the skin.dat from Andrews and Herzberg (1985) book's year,
%             melanoma and sunspot number series.
%
% Syntax:
%   u = SUNMELADAT(selection);
%
% Input Argument:
%   selection:  four element row vector with 0 or 1
%               [year male_melanoma total_melanoma sunspot_number]
%
% Output Argument:
%   u:          matrix with Andrews and Herzberg data columns
%
% Examples:
%   sunmeladat([1 1 1 1])  % or 
%   sunmeladat()           % return all series
%
%   sunmeladat([0 0 1 1])  % return total melanoma and sunspot series
%
%   sunmeladat([4 3]) % return sunspot and melanoma series in this order
%
%   sunmeladat(1)     % or
%   sunmeladat([1])   % yield year column
%
% Reference:
%   Data sample borrowed from 
%      Andrews DF, Herzberg AM. (1985) Data: A Collection of Problems from 
%      Many Fields for the Student and Research Worker. Springer, New York.
%      <https://doi.org/10.1007/978-1-4612-5098-2>
%      ISBN: 978-1-4612-9563-1 (Print) 978-1-4612-5098-2 (Online) 
%          
% Description
%  "This data is from skin.dat Data File in Andrews and Herzberg (1985).
%   The aetiology of melanoma is complex and may include the influences 
%   of trauma, heredity and hormonal activity. In particular, exposure 
%   to solar radiation may be involved in the pathogenesis of melanoma. 
%   Melanoma is more common in fair-skinned individuals and most frequent 
%   in skin sites exposed to the sun. In white populations melanoma is 
%   more common in areas closer to the equator where the intensity of 
%   solar radiation is higher. Data from various parts of the world suggest 
%   that the incidence of melanoma is increasing. The data below, giving 
%   age-adjusted melanoma incidence, are from the Connecticut Tumor Registry 
%   from 1936-1972. Connecticut has the longest record of state 
%   population-based cancer statistics in the United States of America. 
%   The data also includes the sunspot relative number. Houghton, Munster 
%   and Viola (1978) have shown that the age-adjusted incidence rate for 
%   malignant melanoma in the state of Connecticut has risen since 1935 and 
%   that superimposed on the rise are 3-5 year periods in which the rise 
%   in the rate of incidence is excessive. These periods have a cycle of 
%   8-11 years and follow times of maximum sunspot activity. The relationship 
%   between solar cycles and melanoma supports the hypothesis that melanoma 
%   is related to sun exposure and provides evidence that solar radiation 
%   may trigger the development of clinically apparent melanoma. The columns 
%   are the year, male incidence, total incidence, and sunspot relative 
%   index. The incidence are rates per 100,000." (Andrews & Herzberg, 1985)
% 
%   Below is the contents of the file skin.dat:
%   Year % 1936-1972
%   ANNUAL MALE MELANOMA INCIDENCE(AGE-ADJUSTED PER 10**5) CONNECTICUT
%   ANNUAL TOTAL MELANOMA INCIDENCE(AGE-ADJUSTED PER 10**5) CONNECTICUT
%   ANNUAL SUNSPOT RELATIVE NUMBER 


   y=[1936   1.0   0.9    40;
      1937   0.8   0.8   115;
      1938   0.8   0.8   100;
      1939   1.4   1.3    80;
      1940   1.2   1.4    60;
      1941   1.0   1.2    40;
      1942   1.5   1.7    23;
      1943   1.9   1.8    10;
      1944   1.5   1.6    10;
      1945   1.5   1.5    25;
      1946   1.5   1.5    75;
      1947   1.6   2.0   145;
      1948   1.8   2.5   130;
      1949   2.8   2.7   130;
      1950   2.5   2.9    80;
      1951   2.5   2.5    65;
      1952   2.4   3.1    20;
      1953   2.1   2.4    10;
      1954   1.9   2.2     5;
      1955   2.4   2.9    10;
      1956   2.4   2.5    60;
      1957   2.6   2.6   190;
      1958   2.6   3.2   180;
      1959   4.4   3.8   175;
      1960   4.2   4.2   120;
      1961   3.8   3.9    50;
      1962   3.4   3.7    35;
      1963   3.6   3.3    20;
      1964   4.1   3.7    10;
      1965   3.7   3.9    15;
      1966   4.2   4.1    30;
      1967   4.1   3.8    60;
      1968   4.1   4.7   105;
      1969   4.0   4.4   105;
      1970   5.2   4.8   105;
      1971   5.3   4.8    80;
      1972   5.3   4.8    65];

   choice = [];

   if nargin == 0
      choice = [1 2 3 4];
   else
      if max(selection) == 1,
         for i = 1:max(size(selection)),
            if selection(i)
               choice = [choice i];
            end;
         end;
      elseif max(selection) < 1,
         error('Choose at least one column of Andrews and Herzberg data.')
      else
         choice = selection;
      end
   end;

u = y(:,choice);

% These data correspond to Table 32.1 of
% http://www.isds.duke.edu/courses/Spring01/sta114/data/andrews.html
