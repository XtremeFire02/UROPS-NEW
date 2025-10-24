CONFIG -= qt thread import_plugins app_bundle lib_bundle
CONFIG += release optimize_full c++2b
CONFIG += console

TEMPLATE = app

TARGET = shipR

SOURCES += \
    csvFn.cpp \
    ship.cpp \
    shipBioFoul.cpp \
    shipHm.cpp \
    shipResist.cpp \
    shipTable.cpp \
    shipV.cpp \
    prop.cpp \
    wave.cpp \
    waveRlp.cpp \
    waveRcth.cpp \
    wind.cpp \

HEADERS += \
    csvFn.h \
    ship.h \
    shipV.h \
    prop.h \
    wave.h \
    wind.h \
    rotorF.h \

DESTDIR = $$PWD/_bin

unix {
    QMAKE_LFLAGS += -static-libstdc++

    QMAKE_CXXFLAGS += -ffast-math #for pow(x, 2), similarly for clang

    QMAKE_CXXFLAGS += -std=c++2b
    QMAKE_CXXFLAGS += -fopenmp
    QMAKE_LFLAGS += -fopenmp
}

INCLUDEPATH += \
    ../ \

QMAKE_TARGET.arch = $$QMAKE_HOST.arch
DEFINES *= NDEBUG
