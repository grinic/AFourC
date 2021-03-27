#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 10:57:50 2019

@author: ngritti
"""
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QApplication, QComboBox, QVBoxLayout,
        QGridLayout, QGroupBox, QLabel, QPushButton, QDoubleSpinBox,
        QFileDialog, QWidget, QMessageBox, QCheckBox,
        QSpinBox, QListView, QTreeView, QFileSystemModel,
        QAbstractItemView, QDialog)
# from fbs_runtime.application_context.PyQt5 import ApplicationContext
import sys, warnings, os
import viewer4C
import overview4C
warnings.filterwarnings("ignore")


class WindowViewer(QWidget):
    def __init__(self, files, parent=None):
        super(WindowViewer, self).__init__(parent)

        self.files = files

        self.chromosome = QComboBox()
        for i in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y','M']:
            self.chromosome.addItem('chr'+str(i))
        self.chromosome.setCurrentIndex(0)
        
        self.XaxLimStart = QSpinBox()
        self.XaxLimStart.setMaximum(2147483647)
        self.XaxLimStart.setValue(0)
        self.XaxLimStart.resize(50,10)

        self.XaxLimStop = QSpinBox()
        self.XaxLimStop.setMaximum(2147483647)
        self.XaxLimStop.setValue(2147483647)
        self.XaxLimStart.resize(50,10)

        self.closeContact = QSpinBox()
        self.closeContact.setMaximum(2147483647)
        self.closeContact.setValue(35000)
        self.closeContact.resize(50,10)

        self.pValue = QDoubleSpinBox()
        self.pValue.setDecimals(8)
        self.pValue.setValue(0.0005)
        self.pValue.resize(50,10)

        self.windowSize = QSpinBox()
        self.windowSize.setMaximum(10000)
        self.windowSize.setValue(30)
        self.windowSize.resize(50,10)

        self.filesLbl = []
        self.vp_start = []
        self.vp_stop = []
        self.background = []
        self.peakProminence = []
        self.peakWidth = []
        for i in files:
            fileCheck = QCheckBox(os.path.split(i)[-1])
            self.filesLbl.append(fileCheck)
            vp = QSpinBox()
            vp.setMaximum(2147483647)
            vp.setValue(2147483647/2)
            vp.resize(50,10)
            self.vp_start.append(vp)
            vp = QSpinBox()
            vp.setMaximum(2147483647)
            vp.setValue(2147483647/2)
            vp.resize(50,10)
            self.vp_stop.append(vp)
            spinbox = QSpinBox()
            spinbox.setMaximum(10000)
            spinbox.setValue(200)
            spinbox.resize(50,3)
            self.background.append(spinbox)
            spinbox = QSpinBox()
            spinbox.setMaximum(2147483647)
            spinbox.setValue(200)
            spinbox.resize(50,3)
            self.peakProminence.append(spinbox)
            spinbox = QSpinBox()
            spinbox.setMaximum(2147483647)
            spinbox.setValue(16)
            spinbox.resize(50,3)
            self.peakWidth.append(spinbox)

        visualButton = QPushButton('Visualize 4C data')
        visualButton.clicked.connect(self.run_visualization)

        self.strengthLowerLimit = QDoubleSpinBox()
        self.strengthLowerLimit.setMaximum(1)
        self.strengthLowerLimit.setValue(0.75)
        self.strengthLowerLimit.resize(50,10)

        visualMutualButton = QPushButton('Visualize mutual 4C data')
        visualMutualButton.clicked.connect(self.run_mutual_visualization)

        self.verboseButton = QCheckBox('Verbose plotting')
        self.colorbarButton = QCheckBox('Plot colorbar')

        # create the file group
        group = QGroupBox("")
        layout_g = QGridLayout()
        layout_g.addWidget(QLabel('4C files:'),             0,0,1,1)
        layout_g.addWidget(QLabel('Viewpoint start:'),      0,1,1,1)
        layout_g.addWidget(QLabel('Viewpoint stop:'),       0,2,1,1)
        layout_g.addWidget(QLabel('Background:'),           0,3,1,1)
        layout_g.addWidget(QLabel('Peak prominence:'),      0,4,1,1)
        layout_g.addWidget(QLabel('Peak width:'),           0,5,1,1)
        i=1
        for f, vstart, vstop, bckg, pp, pw in zip(self.filesLbl,self.vp_start,self.vp_stop, self.background, self.peakProminence, self.peakWidth):
            layout_g.addWidget(f,i,         0,1,1)
            layout_g.addWidget(vstart,i,    1,1,1)
            layout_g.addWidget(vstop,i,     2,1,1)
            layout_g.addWidget(bckg,i,      3,1,1)
            layout_g.addWidget(pp,i,        4,1,1)
            layout_g.addWidget(pw,i,        5,1,1)
            i+=1

        automaticBackground = QPushButton('Estimate parameters')
        automaticBackground.clicked.connect(self.estimate_parameters)
        layout_g.addWidget(automaticBackground,            i,0,1,6)
        group.setLayout(layout_g)
        ###

        layout = QGridLayout(self)
        layout.addWidget(QLabel('Chromosome'),                      0,0,1,1)
        layout.addWidget(self.chromosome,                           0,1,1,1)
        layout.addWidget(QLabel('Start'),                           1,0,1,1)
        layout.addWidget(self.XaxLimStart,                          1,1,1,1)
        layout.addWidget(QLabel('End'),                             2,0,1,1)
        layout.addWidget(self.XaxLimStop,                           2,1,1,1)
        layout.addWidget(QLabel('Close contact region (bp)'),       3,0,1,1)
        layout.addWidget(self.closeContact,                         3,1,1,1)
        layout.addWidget(QLabel('p-value'),                         4,0,1,1)
        layout.addWidget(self.pValue,                               4,1,1,1)
        layout.addWidget(QLabel('Running window size (fragments)'), 5,0,1,1)
        layout.addWidget(self.windowSize,                           5,1,1,1)
        layout.addWidget(self.verboseButton,                        6,0,1,1)
        layout.addWidget(self.colorbarButton,                       7,0,1,1)
        layout.addWidget(group,                                     8,0,1,2)
        layout.addWidget(visualButton,                              9,0,1,2)
        layout.addWidget(QLabel('Lower limit strength detection'),  10,0,1,1)
        layout.addWidget(self.strengthLowerLimit,                   10,1,1,1)
        layout.addWidget(visualMutualButton,                        11,0,1,2)
        self.setLayout(layout)

        self.setWindowTitle("AFourC - Analysis")

    def estimate_parameters(self):
        files = self.files
        chromosome = self.chromosome.currentText()

        for i, file in enumerate(files):
            background, prominence = viewer4C.estimate_parameters(file, chromosome)
            self.background[i].setValue(background)
            self.peakProminence[i].setValue(prominence)

    def run_visualization(self):
        print('running visualization...')
        N_files = len(self.files)
        chromosome = self.chromosome.currentText()
        _lims = (self.XaxLimStart.value(),self.XaxLimStop.value())

        files = [self.files[i] for i in range(N_files) if self.filesLbl[i].isChecked()]
        viewpoints = [[self.vp_start[i].value(),self.vp_stop[i].value()] for i in range(N_files) if self.filesLbl[i].isChecked()]
        backgrounds = [self.background[i].value() for i in range(N_files) if self.filesLbl[i].isChecked()]
        peakproms = [self.peakProminence[i].value() for i in range(N_files) if self.filesLbl[i].isChecked()]
        peakwidths = [self.peakWidth[i].value() for i in range(N_files) if self.filesLbl[i].isChecked()]

        result = viewer4C.plot_4C_interactions(
                        files, 
                        chromosome, 
                        viewpoints,
                        _lims, 
                        close_contact_region=self.closeContact.value(),
                        colors=None, 
                        colors_profile=None,
                        colorbar=self.colorbarButton.isChecked(), 
                        seed=1,
                        windowSize=self.windowSize.value(),
                        backgrounds=backgrounds,
                        peaks_prominence=peakproms,
                        peaks_width=peakwidths,
                        pvalue=self.pValue.value(),
                        verbose=self.verboseButton.isChecked()
                        )
        if result == False:
            QMessageBox.warning(self,'Warning!','I couldn\'t find any data in this Chromosome or region!')

    def run_mutual_visualization(self):
        print('running mutual visualization...')
        N_files = len(self.files)
        if N_files < 2:
            QMessageBox.warning(self,'Warning!','You need to specify at least two files to copute mutual interaction!')
            return

        files = self.files
        chromosome = self.chromosome.currentText()
        viewpoints = [[self.vp_start[i].value(),self.vp_stop[i].value()] for i in range(N_files)]
        _lims = (self.XaxLimStart.value(),self.XaxLimStop.value())
        close_contact_region = self.closeContact.value()
        backgrounds = [self.background[i].value() for i in range(N_files) if self.filesLbl[i].isChecked()]
        peakproms = [self.peakProminence[i].value() for i in range(N_files) if self.filesLbl[i].isChecked()]
        peakwidths = [self.peakWidth[i].value() for i in range(N_files) if self.filesLbl[i].isChecked()]

        result = overview4C.plot_4C_overview(
                    files, 
                    chromosome, 
                    viewpoints, 
                    _lims, 
                    close_contact_region=close_contact_region,
                    strength_lower_limit=self.strengthLowerLimit.value(),
                    colors=None, 
                    colors_profile=None,
                    colorbar=self.colorbarButton.isChecked(), 
                    seed=1,
                    windowSize=self.windowSize.value(),
                    backgrounds=backgrounds,
                    peaks_prominence=peakproms,
                    peaks_width=peakwidths,
                    pvalue=self.pValue.value(),
                    verbose=self.verboseButton.isChecked()
                    )

        if result == False:
            QMessageBox.warning(self,'Warning!','Something went wrong!')
                 
class FileDialog(QFileDialog):
    def __init__(self, *args):
        QFileDialog.__init__(self, *args)
        self.setOption(self.DontUseNativeDialog, True)
        # self.setFileMode(self.DirectoryOnly)

        for view in self.findChildren((QListView, QTreeView)):
            if isinstance(view.model(), QFileSystemModel):
                view.setSelectionMode(QAbstractItemView.ExtendedSelection)

class app4C(QWidget):
    def __init__(self, parent=None):
        super(app4C, self).__init__(parent)

        selectBtn = QPushButton('Select 4C input files')
        selectBtn.clicked.connect(self.selectFiles)

        layout = QVBoxLayout(self)
        layout.addWidget(selectBtn)
        self.setLayout(layout)

        self.setWindowTitle("AFourC - Starting window")
        # self.setSize(300, 85)

    def selectFiles(self):
        dialog = FileDialog()
        if dialog.exec_() == QDialog.Accepted:
            files = dialog.selectedFiles()
        else:
            return

        # inds = self.tree.selectionModel().selectedIndexes()
        # files = []
        # for i in inds:
        #     if i.column() == 0:
        #         files.append(os.path.join(str(self.directory().absolutePath()),str(i.data().toString())))
        # self.selectedFiles = files
        # self.hide()

        self.w = WindowViewer(files)
        self.w.resize(200,60)
        self.w.show()
    # def selectFiles(self):
    #     file_name = QFileDialog()
    #     file_name.setFileMode(QFileDialog.ExistingFiles)
    #     files, _ = file_name.getOpenFileNames(self, "Open files", "C\\Desktop")
    #     print(files)

    #     self.w = WindowViewer(files)
    #     self.w.show()
        
if __name__ == '__main__':
    appctxt = QApplication(sys.argv)       # 1. Instantiate ApplicationContext
    window = app4C()
    window.resize(400, 80)
    window.show()
    exit_code = appctxt.exec_()      # 2. Invoke appctxt.app.exec_()
    sys.exit(exit_code)
