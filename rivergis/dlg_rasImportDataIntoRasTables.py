# -*- coding: utf-8 -*-

import psycopg2
import psycopg2.extras

from PyQt4.QtCore import *
from PyQt4.QtGui import *
from qgis.core import *
from qgis.gui import *
from qgis.utils import *
from ui.ui_importDataIntoRasTables import *


class DlgImportDataIntoRasTables(QDialog):
  def __init__(self, rgis):
    QDialog.__init__(self)
    self.ui = Ui_importDataIntoRasTables()
    self.ui.setupUi(self)
    self.rgis = rgis
    self.ui.buttonBox.accepted.connect(self.accept)
    self.ui.buttonBox.rejected.connect(self.rejectDlg)
    self.ui.cboBanks.currentIndexChanged.connect(self.cboBanksLayerChanged)
    self.ui.cboFlowPaths.currentIndexChanged.connect(self.cboFlowpathsLayerChanged)
    self.ui.cboIneffective.currentIndexChanged.connect(self.cboIneffectiveLayerChanged)
    self.ui.cboObstructions.currentIndexChanged.connect(self.cboObstructionsLayerChanged)
    # QObject.connect(self.ui.helpButton, SIGNAL("clicked()"), self.displayHelp)
    self.populateCbos()
    self.banksLayer = None
    self.flowpathsLayer = None
    self.ineffectiveLayer = None
    self.obstructionsLayer = None


  def accept(self):
    QApplication.setOverrideCursor(Qt.WaitCursor)

    if not self.ui.cboRivers.currentText() == '':
        # TODO: import warstwy rzek do PG
        pass

    if not self.ui.cboXsecs.currentText() == '':
        # TODO: import warstwy przekrojow do PG
        pass

    if self.banksLayer:
      # TODO: import warstwy banksLayer do PG

      if not self.ui.cboBanklineType.currentText() == '': # Bank type (side) specified
        # TODO: kod uwzgledniajacy atrybut typu (strony) brzegu
        pass

    if self.flowpathsLayer:
      # TODO: import warstwy linii przeplywow do PG

      if not self.ui.cboFlowpathType.currentText() == '': # Bank type (side) specified
        # TODO: kod uwzgledniajacy atrybut typu linii
        pass

    if not self.ui.cboLevees.currentText() == '':
        # TODO: import warstwy walow do PG
        pass

    if self.ineffectiveLayer:
      # TODO: import warstwy pol jalowego przeplywu do PG

      if not self.ui.cboFlowpathType.currentText() == '': # Bank type (side) specified
        # TODO: kod uwzgledniajacy atrybut wysokosci pola jalowego
        pass

    if self.obstructionsLayer:
      # TODO: import warstwy przeszkod do PG

      if not self.ui.cboObstructionsElev.currentText() == '': # Bank type (side) specified
        # TODO: kod uwzgledniajacy atrybut wysokosci przeszkody
        pass

    self.rgis.addInfo("  Import completed.")
    QApplication.setOverrideCursor(Qt.ArrowCursor)
    QDialog.accept(self)




  def displayHelp(self):
    pass

  def populateCbos(self):
    self.ui.cboRivers.clear()
    self.ui.cboXs.clear()
    self.ui.cboLevees.clear()
    self.ui.cboBanks.clear()
    self.ui.cboIneffective.clear()
    self.ui.cboObstructions.clear()
    self.ui.cboFlowPaths.clear()
    self.ui.cboRivers.addItem("")
    self.ui.cboXs.addItem("")
    self.ui.cboLevees.addItem("")
    self.ui.cboBanks.addItem("")
    self.ui.cboIneffective.addItem("")
    self.ui.cboObstructions.addItem("")
    self.ui.cboFlowPaths.addItem("")
    for layerId, layer in self.rgis.mapRegistry.mapLayers().iteritems():
      if layer.type() == 0 and layer.geometryType() == 0: # vector and points
        pass
      if layer.type() == 0 and layer.geometryType() == 1: # vector and polylines
        self.ui.cboRivers.addItem(layer.name(), layerId)
        self.ui.cboXs.addItem(layer.name(), layerId)
        self.ui.cboLevees.addItem(layer.name(), layerId)
        self.ui.cboBanks.addItem(layer.name(), layerId)
        self.ui.cboFlowPaths.addItem(layer.name(), layerId)
      if layer.type() == 0 and layer.geometryType() == 2: # vector and polygons
        self.ui.cboIneffective.addItem(layer.name(), layerId)
        self.ui.cboObstructions.addItem(layer.name(), layerId)
      if layer.type() == 1: # it's a raster
        pass


  def updateBanksAttrs(self):
    if self.banksLayer:
      if self.banksLayer.featureCount():
        self.ui.cboBanklineType.clear()
        self.ui.cboBanklineType.addItem("")
        attrs = self.banksLayer.pendingFields()
        typeIdx = 0
        for i, attr in enumerate(attrs):
          self.ui.cboBanklineType.addItem(attr.name())
          # check if there is 'cellsize' attr
          if attr.name() == 'type' or attr.name() == 'typ':
            typeIdx = i + 1
        self.ui.cboBanklineType.setCurrentIndex(typeIdx)


  def cboBanksLayerChanged(self):
    curInd = self.ui.cboBanks.currentIndex()
    lid = self.ui.cboBanks.itemData(curInd)
    self.banksLayer = self.rgis.mapRegistry.mapLayer(lid)
    self.updateBanksAttrs()


  def updateFlowpathsAttrs(self):
    if self.flowpathsLayer:
      if self.flowpathsLayer.featureCount():
        self.ui.cboFlowpathType.clear()
        self.ui.cboFlowpathType.addItem("")
        attrs = self.flowpathsLayer.pendingFields()
        typeIdx = 0
        for i, attr in enumerate(attrs):
          self.ui.cboFlowpathType.addItem(attr.name())
          if attr.name().lower() == 'type' or attr.name().lower() == 'typ':
            typeIdx = i + 1
        self.ui.cboFlowpathType.setCurrentIndex(typeIdx)


  def cboFlowpathsLayerChanged(self):
    curInd = self.ui.cboFlowPaths.currentIndex()
    lid = self.ui.cboFlowPaths.itemData(curInd)
    self.flowpathsLayer = self.rgis.mapRegistry.mapLayer(lid)
    self.updateFlowpathsAttrs()


  def updateIneffectiveAttrs(self):
    if self.ineffectiveLayer:
      if self.ineffectiveLayer.featureCount():
        self.ui.cboIneffElev.clear()
        self.ui.cboIneffElev.addItem("")
        attrs = self.ineffectiveLayer.pendingFields()
        typeIdx = 0
        for i, attr in enumerate(attrs):
          self.ui.cboIneffElev.addItem(attr.name())
          if attr.name().lower() == 'elevation' or attr.name().lower() == 'elev':
            typeIdx = i + 1
        self.ui.cboIneffElev.setCurrentIndex(typeIdx)


  def cboIneffectiveLayerChanged(self):
    curInd = self.ui.cboIneffective.currentIndex()
    lid = self.ui.cboIneffective.itemData(curInd)
    self.ineffectiveLayer = self.rgis.mapRegistry.mapLayer(lid)
    self.updateIneffectiveAttrs()


  def updateObstructionsAttrs(self):
    if self.obstructionsLayer:
      if self.obstructionsLayer.featureCount():
        self.ui.cboObstructionsElev.clear()
        self.ui.cboObstructionsElev.addItem("")
        attrs = self.obstructionsLayer.pendingFields()
        typeIdx = 0
        for i, attr in enumerate(attrs):
          self.ui.cboObstructionsElev.addItem(attr.name())
          if attr.name().lower() == 'elevation' or attr.name().lower() == 'elev':
            typeIdx = i + 1
        self.ui.cboObstructionsElev.setCurrentIndex(typeIdx)


  def cboObstructionsLayerChanged(self):
    curInd = self.ui.cboObstructions.currentIndex()
    lid = self.ui.cboObstructions.itemData(curInd)
    self.obstructionsLayer = self.rgis.mapRegistry.mapLayer(lid)
    self.updateObstructionsAttrs()


  def rejectDlg(self):
    self.rgis.addInfo("  Importing data cancelled.")
    self.reject()



