function    cssm( filename, linenum ) 
    % deletes the linenum-th line in the file given by filename
    
    
    str = fileread( filename ); 
    pos = strfind( str, sprintf('\n') );
    str( pos(linenum-1)+1 : pos(linenum) ) = [];
    fid = fopen( strcat('contig_',filename), 'w' );
    fprintf( fid, '%s', str );
    fclose( fid );
end
