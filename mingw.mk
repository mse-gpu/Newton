# Version 	: 0.0.1
# Date		: 26.01.2012
# Author 	: Cedric.Bilat@he-arc.ch
#
# Contrainte
#
#	Dans les définitions de variables ci-dessous, ne laisser jamais d'espace à la fin!!!
#

#########
#public #
#########

SRC_AUX+=../PRODUCTION/commonExt/include
SRC_AUX+=../PRODUCTION/commonExt/windows

#Option de compilation 
#Elles seront ensuite surchargées automatiquement (par exemple en fonction du type de la target)
#Compilation: (minimal conseillé : vide) 
#CXXFLAGS:= -W -Wall 
#CXXFLAGS+= -mtune=core2 
CXXFLAGS+= -mtune=corei7
CXXFLAGS+= -O3 -frerun-loop-opt -fgcse -fomit-frame-pointer
CXXFLAGS+= -fopenmp

#Option de compilation Static Lib 
ARFLAGS:=
 
#Option de linkage
#Elles seront ensuite surchargées automatiquement (par exemple en fonction du type de la target)
#Version minimale : Ne rien metter laisser vierge!
LDFLAGS:=

#Injection de variable dans le code (same as #define XXX YYY)
CODE_DEFINE_VARIABLES:=

#separateur est espace,sans extension .lib, sans prefixe lib
ADD_LIBRARY_FILES:=

EXCLUDE_LIBRARY_FILES:=

#########
#Private#
#########

COMPILATEUR:=mingw
OS:=Linux
USER=${USERNAME}

ROOT_MAKEFILE_PRIVATE:=../PRODUCTION/MANAGER/makefile/private
-include dataProject.mk
-include $(ROOT_MAKEFILE_PRIVATE)/fileSet.mk
-include $(ROOT_MAKEFILE_PRIVATE)/flags_MinGW.mk
-include $(ROOT_MAKEFILE_PRIVATE)/makeFileCPP.mk
 
#########
#  End  #
#########