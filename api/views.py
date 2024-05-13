from core.models import Category, Molecule
from core.views import DefaultPagination
from django.db.models import Count
from rest_framework import filters, generics, pagination, viewsets, status
from rest_framework.decorators import action
from rest_framework.response import Response
from .filters import MoleculeFilter
from django_filters.rest_framework import DjangoFilterBackend

from core.serializers import CategorySerializer, MoleculeSerializer

import time
import os


class OverviewViewSet(viewsets.ReadOnlyModelViewSet):
    queryset = Molecule.objects.all()
    serializer_class = MoleculeSerializer
    pagination_class = DefaultPagination


class DownloadViewSet(viewsets.ModelViewSet):
    queryset = Molecule.objects.all()
    serializer_class = MoleculeSerializer
    pagination_class = DefaultPagination

    @action(detail=False, methods=["get"], url_path="all")
    def all_classes(self, request):
        queryset = Molecule.objects.all()

    @action(detail=False, methods=["get"], url_path="all")
    def all_classes(self, request):
        categories = Category.objects.all()
        paginated_categories = self.paginate_queryset(
            categories
        )  # Apply pagination to categories
        result = []

        if paginated_categories is not None:
            for category in paginated_categories:
                molecules = Molecule.objects.filter(class_type=category)
                serializer = self.get_serializer(molecules, many=True)
                result.append(
                    {"class_type": category.name, "molecules": serializer.data}
                )

            return self.get_paginated_response(result)  # Return paginated response
        else:
            # Handle the case if pagination is not applicable or error
            return Response({"error": "Pagination error or no categories found."})


# api for searching with multiple query params single molecule data
class SearchViewSet(viewsets.ModelViewSet):
    queryset = Molecule.objects.all()
    serializer_class = MoleculeSerializer  # TODO: this serves search api, but should have its own serializer
    pagination_class = DefaultPagination

    filter_backends = (DjangoFilterBackend,)
    filterset_class = MoleculeFilter


# api for category object
class CategoryViewSet(viewsets.ViewSet):
    def list(self, request):
        categories = Category.objects.all()
        serializer = CategorySerializer(categories, many=True)
        return Response(serializer.data)

    def retrieve(self, request, pk=None):
        molecules = Molecule.objects.filter(class_type__id=pk)
        serialized_molecules = MoleculeSerializer(molecules, many=True)
        return Response(serialized_molecules.data)


# api for statistics
class StatisticsViewSet(viewsets.ViewSet):
    def list(self, request):
        total_molecules = Molecule.objects.count()
        total_categories = Category.objects.count()

        return Response(
            {
                "total_molecules": total_molecules,
                "total_categories": total_categories,
            }
        )


# api for weight distribution stats
class WeightDistributionViewSet(viewsets.ViewSet):
    def weight_distribution(self):
        weight_ranges = ["0-100", "100-200", "200-300", "300-400", "400-500", "500+"]

        weight_bins = {
            "0-100": 0,
            "100-200": 0,
            "200-300": 0,
            "300-400": 0,
            "400-500": 0,
            "500+": 0,
        }

        # Adjust weights and classify into bins
        for molecule in Molecule.objects.all():
            weight = molecule.molecular_weight
            if weight is not None:
                if 0 < weight < 100:
                    weight_bins["0-100"] += 1
                elif weight < 200:
                    weight_bins["100-200"] += 1
                elif weight < 300:
                    weight_bins["200-300"] += 1
                elif weight < 400:
                    weight_bins["300-400"] += 1
                elif weight < 500:
                    weight_bins["400-500"] += 1
                else:
                    weight_bins["500+"] += 1

        # Prepare data for the chart
        data = {
            "labels": weight_ranges,
            "values": [weight_bins[range_] for range_ in weight_ranges],
        }

        return data

    def list(self, request):
        return Response(self.weight_distribution())


class ChiralityDistributionViewSet(viewsets.ViewSet):
    def get_chirality_distribution(self):
        distribution_data = list(
            Molecule.objects.values("chirality__name")
            .annotate(count=Count("id"))
            .order_by()
        )

        return distribution_data

    def list(self, request):
        return Response(self.get_chirality_distribution())


class CategoryDistributionViewSet(viewsets.ViewSet):
    def get_category_distribution(self):
        distribution_data = list(
            Molecule.objects.values("category__name")
            .annotate(count=Count("id"))
            .order_by()
        )
        return distribution_data

    def list(self, request):
        distribution_data = self.get_category_distribution()
        return Response(distribution_data)


class SDFUploaderViewSet(viewsets.ViewSet):
    def create(self, request):
        uploaded_files = request.FILES.getlist("files")

        if not uploaded_files:
            return Response(
                {"error": "No files uploaded."}, status=status.HTTP_400_BAD_REQUEST
            )

        random_folder = str(time.time()).replace(".", "-")
        os.makedirs(f"cluster/{random_folder}", exist_ok=True)

        # Process and save the uploaded files
        saved_files = []
        for file in uploaded_files:
            # Save the file to a desired location

            file_path = f"cluster/{random_folder}/{file.name}"
            with open(file_path, "wb") as f:
                f.write(file.read())
            saved_files.append(file.name)

        # Log the received files
        print(f"Received files: {saved_files}")

        return Response(
            {"message": "Files uploaded successfully", "files": saved_files},
            status=status.HTTP_200_OK,
        )
