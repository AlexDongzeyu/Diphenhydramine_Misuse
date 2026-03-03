March 02, 2026

README: RxNorm 03/02/2026 Current Prescribable Content 
===================================================

-----------------------------------------------------------------
This current prescribable content release contains data that is consistent with the
2020AA version of the UMLS.
-----------------------------------------------------------------
This release has been created from the March 02, 2026 Full
Release version of RxNorm.

For full details, please refer to the RxNorm documentation at
https://www.nlm.nih.gov/research/umls/rxnorm/docs/prescribe.html.

This release contains database control files and SQL
commands for use in the automation of the loading process of
these files into an Oracle or MySQL database.

RxNorm release data files are available by download from
the NLM download server at:

        https://www.nlm.nih.gov/research/umls/rxnorm/docs/rxnormfiles.html

This link will take you to a page for downloading the latest files (no login required):
RxNorm_full_prescribe_03022026.zip

Once downloaded, it must be unzipped in order to access the files.

HARDWARE AND SOFTWARE RECOMMENDATIONS
-------------------------------------
- Supported operating systems:
        Windows: 7
        Linux
        Solaris: Solaris 10

- Hardware Requirements

  - A MINIMUM 560 MB of free hard disk space (To accomodate ZIP files and
        unzipped contents).

CONTENTS OF THE ZIP FILE
-------------------------

The ZIP formatted file contains the
following 16 files and 4 directories:

Readme_Full_Prescribe_03022026.txt  1031                 bytes

rrf directory:

RXNCONSO.RRF                        30,367,709           bytes
RXNREL.RRF                          193,258,180          bytes
RXNSAT.RRF                          275,372,939          bytes


scripts directory:

        oracle sub-directory:

populate_oracle_rxn_db.bat          699                  bytes
RXNCONSO.ctl                        512                  bytes
rxn_index.sql                       460                  bytes
RxNormDDL.sql                       1,373                bytes
RXNREL.ctl                          471                  bytes
RXNSAT.ctl                          378                  bytes

         mysql sub-directory:

Indexes_mysql_rxn.sql               463                  bytes
Load_scripts_mysql_rxn_unix.sql     1,469                bytes
Load_scripts_mysql_rxn_win.sql      1,468                bytes
Populate_mysql_rxn.bat              777                  bytes
populate_mysql_rxn.sh               1,609                bytes
Table_scripts_mysql_rxn.sql         1,749                bytes
Additional NOTES:
-----------------

- Most RxNorm users will need applications and data management
  systems such as an RDBMS for storage and retrieval.

- The RxNorm release files contain UTF-8 Unicode encoded data.

- Refer to the RxNorm prescribing release documentation at https://www.nlm.nih.gov/research/umls/rxnorm/docs/prescribe.html
  for information on the contents of the RxNorm Current Prescribable Content release files.
