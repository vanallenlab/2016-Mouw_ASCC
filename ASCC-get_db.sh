#!/bin/bashß

# Usage
# bash ASCC-get_db.sh

echo 'Running ASCC-get_db.sh...'

echo 'Download COSMIC v75 from Dropbox' 
# Using curl. If on linux, use wget
curl -L -o cosmic_v75.txt https://www.dropbox.com/sh/jm3jqgzrlgmshmp/cosmic_v75.txt?dl=1

echo 'Downloading ExAC from Dropbox...'
curl -L -o exac_GATK.txt https://www.dropbox.com/sh/hitimba3ra0vdgb/exac_GATK.txt?dl=1

echo 'Moving COSMIC v75 and ExAC to appropriate directory'
mv cosmic_v75.txt storage/reference/cosmic_v75.txt
mv exac_GATK.txt storage/reference/exac_GATK.txt

# References:
# 
# COSMIC: exploring the world's knowledge of somatic mutations in human cancer (Forbes et al. 2014)
# Cosmic(v75)
# http://cancer.sanger.ac.uk/cosmic
# 
# http://exac.broadinstitute.org
# Exome Aggregation Consortium, M. Lek, K. Karczewski, E. Minikel, K. Samocha, E. Banks, T. Fennell, A. O’Donnell-Luria, J. Ware, A. Hill, B. Cummings, T. Tukiainen, D. Birnbaum, J. Kosmicki, L. Duncan, K. Estrada, F. Zhao, J. Zou, E. Pierce-Hoffman, D. Cooper, M. DePristo, R. Do, J. Flannick, M. Fromer, L. Gauthier, J. Goldstein, N. Gupta, D. Howrigan, A. Kiezun, M. Kurki, A. L. Moonshine, P. Natarajan, L. Orozco, G. Peloso, R. Poplin, M. Rivas, V. Ruano-Rubio, D. Ruderfer, K. Shakir, P. Stenson, C. Stevens, B. Thomas, G. Tiao, M. Tusie-Luna, B. Weisburd, H.-H. Won, D. Yu, D. Altshuler, D. Ardissino, M. Boehnke, J. Danesh, E. Roberto, J. Florez, S. Gabriel, G. Getz, C. Hultman, S. Kathiresan, M. Laakso, S. McCarroll, M. McCarthy, D. McGovern, R. McPherson, B. Neale, A. Palotie, S. Purcell, D. Saleheen, J. Scharf, P. Sklar, S. Patrick, J. Tuomilehto, H. Watkins, J. Wilson, M. Daly, D. MacArthur, Analysis of protein-coding genetic variation in 60,706 humans. bioRxiv 10.1101/030338 (2015)
# This file was processed using GATK
# https://www.broadinstitute.org/gatk/
