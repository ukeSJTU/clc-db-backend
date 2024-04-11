from core.models import Category, Molecule
from core.views import OverviewPagination
from django.db.models import Prefetch
from rest_framework import filters, generics, pagination, viewsets
from rest_framework.decorators import action
from rest_framework.response import Response

from .serializers import OverviewSerializer


# cards grid layout with pagination
class OverviewCardViewSet(viewsets.ReadOnlyModelViewSet):
    queryset = Molecule.objects.all()
    serializer_class = OverviewSerializer
    pagination_class = OverviewPagination


# table layout with pagination
class OverviewTableView(viewsets.ReadOnlyModelViewSet):
    queryset = Molecule.objects.all()
    serializer_class = OverviewSerializer
    pagination_class = OverviewPagination
