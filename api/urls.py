from django.urls import include, path
from rest_framework.routers import DefaultRouter

from .views import *

router = DefaultRouter()

router.register(r"overview/card", OverviewCardViewSet, basename="overview-card")
router.register(r"overview/table", OverviewTableView, basename="overview-table")


urlpatterns = [
    path("", include(router.urls)),
]
