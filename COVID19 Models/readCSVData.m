function [dates_extended, newcases_extended, totalcases_extended] = readCSVData(CSVFile, CountryNames)
%% reads the data in a table
% if the format of data.csv changes, you need to edit this line:
T = readtable(CSVFile,'Format','%{dd/MM/yyyy}D%f%f%f%f%f%q%q%q%f%q', 'Delimiter', ',', 'HeaderLines', 0, 'ReadVariableNames', true);

T_slice = T((strcmp(T.countriesAndTerritories, CountryNames)),:);
[dates, dates_sorted_ind] = sort(datenum(T_slice.year, T_slice.month, T_slice.day));
newcases = T_slice.cases(dates_sorted_ind);

dates_extended = dates(1):dates(end);
newcases_extended = zeros(size(dates_extended));
newcases_extended(arrayfun(@(x)(any(dates==x)), dates_extended))=newcases;
totalcases_extended = cumsum(newcases_extended);

dates_extended = dates_extended(totalcases_extended~=0)';
newcases_extended = newcases_extended(totalcases_extended~=0)';
totalcases_extended = totalcases_extended(totalcases_extended~=0)';

end









