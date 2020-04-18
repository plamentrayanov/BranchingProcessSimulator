function [dates, newcases, totalcases] = readCSVData(CSVFile, CountryNames)
%% reads the data in a table
T = readtable(CSVFile,'Format','%{dd/MM/yyyy}D%f%f%f%f%f%q%q%q%f', 'Delimiter', ',', 'HeaderLines', 0, 'ReadVariableNames', true);
if iscellstr(CountryNames) && length(CountryNames)>1
    for i=1:length(CountryNames)
        T_slice = T((strcmp(T.countriesAndTerritories, CountryNames{i})),:);
        dates{i} = datenum(T_slice.year, T_slice.month, T_slice.day);
        newcases{i} = T_slice.cases;
    end
else
    T_slice = T((strcmp(T.countriesAndTerritories, CountryNames)),:);
    [dates, dates_sorted_ind] = sort(datenum(T_slice.year, T_slice.month, T_slice.day));
    newcases = T_slice.cases(dates_sorted_ind);
    totalcases = cumsum(newcases);
end
end









