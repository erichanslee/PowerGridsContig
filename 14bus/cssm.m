function    cssm( filename, linenum ) 
    %%
    str = fileread( filename ); 
    pos = strfind( str, sprintf('\n') );
    str( pos(linenum-1)+1 : pos(linenum) ) = [];
    %%
    fid = fopen( strcat('contig_',filename), 'w' );
    fprintf( fid, '%s', str );
    fclose( fid );
end
