function TAM_saveCoordinates(filename, x, y)
    % Create a table with x and y
    T = table(x', y', 'VariableNames', {'x', 'y'});
    
    % Check if the file exists
    if isfile(filename)
        % If the file exists, open it in append mode
        %delete(filename); 
        fid = fopen(filename, 'a');
        
        % Write the new data
        fprintf(fid, '%f,%f\n', T{:,:});
    else
        % If the file does not exist, write the table to a new CSV file
        writetable(T, filename);
    end
    
    % Close the file
    if exist('fid', 'var')
        fclose(fid);
    end
end