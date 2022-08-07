QT -= gui

CONFIG += c++11 console
CONFIG -= app_bundle

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

#QMAKE_CXXFLAGS += -flto -W

SOURCES += \
        BarElement_Linear_1D.cpp \
        BarElement_NonLinear_1D.cpp \
        BeamElement_NonLinear_EulerBernoulli_2D.cpp \
        GEBT.cpp \
        ReadInpFile.cpp \
        VAMBeamElement.cpp \
        main.cpp

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

win32:INCLUDEPATH += "G:\MTech_Aerospace\MTech_Res_Thesis\Cpp\TextileComposite\Eigen"

HEADERS += \
    Variables.h
