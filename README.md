# mysql 操作

```bash
brew install mysql
brew update
brew serices start mysql
==> Successfully started `mysql` (label: homebrew.mxcl.mysql)
mysql
```

## To load data

```bash
python manage.py makemigrations
python manage.py migrate
python manage.py load_molecules -p /ChemNexus/data/merged.csv
```
