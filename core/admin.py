from django.contrib import admin

from .models import Category, Molecule, Smile

# Register your models here.
admin.site.register(Molecule)
admin.site.register(Category)
admin.site.register(Smile)
