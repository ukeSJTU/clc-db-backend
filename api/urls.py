from django.urls import include, path
from rest_framework.routers import DefaultRouter

from .views import *

router = DefaultRouter()

router.register(r"overview/card", OverviewViewSet, basename="overview-card")
router.register(r"overview/table", OverviewViewSet, basename="overview-table")

router.register(r"download/molecules", DownloadViewSet, basename="molecule")
router.register(r"download/classes", DownloadViewSet, basename="class")

router.register(r"search/molecules", SearchViewSet, basename="search")

router.register(r"statistics", StatisticsViewSet, basename="statistics")
router.register(r"stats/weights", WeightDistributionViewSet, basename="weight-stats")
router.register(r"stats/smiles", SmileTypeDistributionViewSet, basename="smiles-stats")
router.register(r"stats/classes", ClassTypeDistributionViewSet, basename="class-stats")
router.register(r"categories", CategoryViewSet, basename="category")

urlpatterns = [
    path("", include(router.urls)),
]
