git checkout doc
cp -r build/html/* ..
git add ../*.html
git add ../*.js
git add ../_images/*
git .._static/*
git add ../objects.inv
git commit -m "Updating doc."
git push
git checkout master
