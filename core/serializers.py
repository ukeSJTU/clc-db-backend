from rest_framework import serializers, viewsets

from .models import Category, Chirality


class CategorySerializer(serializers.ModelSerializer):
    class Meta:
        model = Category
        fields = ["name"]


class ChiralitySerializer(serializers.ModelSerializer):
    class Meta:
        model = Chirality
        fields = ["name"]
