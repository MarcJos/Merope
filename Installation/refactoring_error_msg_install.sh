# MJ263790
# Tool for reading more easily the warnings of the compiler, removing some useless data
# Expects : that the output of the compilation is written in erreurs.txt

sed -i '/^--/d' erreurs.txt
sed -i '/voro/,+2 d' erreurs.txt
sed -i '/completeMaterialPropertiesList/d' erreurs.txt
sed -i '/^\[/d' erreurs.txt
sed -i '/Eigen/d' erreurs.txt
sed -i '/LTRANS/d' erreurs.txt
