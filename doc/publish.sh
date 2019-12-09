git checkout gh-pages
cp -r build/* ..
git add ../*.html
git add ../*.js
git add ../_images/*
git add .._static/*
git add ../objects.inv
git commit -m "Updating doc."
git push
git checkout master
