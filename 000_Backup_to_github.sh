# init with 'main' as the default branch
git init -b main

# (optional, if you haven't set these before)
git config --global user.name "Rongjun Huang"
git config --global user.email "u6569836@anu.edu.au"

# create a .gitignore that ignores everything except the types you want
cat > .gitignore << 'EOF'
# Ignore everything by default
*

# But allow directories (so recursion works)
!*/

# Allowlisted file types
!*.ipynb
!*.py
!*.sh
!*.json
!*.log
!*.txt
!*.md
!*.csv
!*.tsv
!*.r
!*.R
!*.Rmd
!*.rmd
!*.yml
!*.yaml
!*.xml
!*.cfg
!*.ini
!*.conf
!*.json
!*.dat
!*.json
!*.dat
!*.png
!*.pdf
!*.tiff
!*.jpg
!*.jpeg
!*.gif

# Keep these metadata files tracked too
!.gitignore
!README.md
!.gitattributes
# Ignore the ppxf folder
ppxf/
Prepare/
Unihall/
# Ignore large files that exceed GitHub limits
extended/MAUVE_all_R.png
MAUVE/Rongjun_rFMR_Virgo_Mauve_Presentation.pdf
extended/MAUVE_all_R.tiff
EOF

# (optional) small README so the repo isn't empty on GitHub
echo "# ICRAR-MAUVE backup" > README.md

git add .
git commit -m "Initial backup of notebooks, scripts, JSON, and logs"


# Then go to GitHub and create a new empty repo called "ICRAR-MAUVE"

# git remote -v                # confirm origin is git@github.com:Rongjun-ANU/ICRAR-MAUVE.git
# git push -u origin main

# git remote set-url origin https://github.com/Rongjun-ANU/ICRAR-MAUVE.git
# git push -u origin main
