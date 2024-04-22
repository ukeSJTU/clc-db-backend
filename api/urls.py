from django.urls import include, path
from rest_framework.routers import DefaultRouter

from .views import *

router = DefaultRouter()

router.register(r"overview/card", OverviewCardViewSet, basename="overview-card")
router.register(r"overview/table", OverviewTableView, basename="overview-table")

router.register(r"download/molecules", DownloadMoleculesViewSet, basename="molecule")
router.register(r"download/classes", DownloadClassesViewSet, basename="class")

router.register(r"search/molecules", MoleculeViewSet, basename="search")


urlpatterns = [
    path("", include(router.urls)),
]
