Identifying Contingencies in the Power Grid using grid dynamics.   

Requires PSAT to be installed and added to system path:  
http://faraday1.ucd.ie/psat.html  
NOTE: PSAT NOT COMPATIBLE WITH MATLAB R2015 AND ABOVE; MUST DOWNGRADE TO AT LEAST R2014  

Each folder will contain a few scripts:  
-genscript.m will generate the contingency files (right now only focused on line failures)  
-calculateEigs.m calculates the eigenvalues of each contingency's state matrix via inverse  
 iteration + Arnoldi + Sherman Morrison Woodbury  
-sanitycheck.m is just a test script to make sure everything works  
-metadata.mat contains metadata for scripts to use  
