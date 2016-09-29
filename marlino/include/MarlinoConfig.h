#ifndef MarlinoConfig_h
#define MarlinoConfig_h 1

// version macros
#define MARLINO_MAJOR_VERSION 1
#define MARLINO_MINOR_VERSION 8
#define MARLINO_PATCH_LEVEL 0

#define MARLINO_VERSION_GE( MAJV , MINV , PLEV )  ( (MARLINO_MAJOR_VERSION > MAJV) || ( (MARLINO_MAJOR_VERSION == MAJV) && (MARLINO_MINOR_VERSION > MINV) ) || ( (MARLINO_MAJOR_VERSION == MAJV) && (MARLINO_MINOR_VERSION == MINV) && (MARLINO_PATCH_LEVEL >= PLEV) ) )

#endif // MarlinConfig_h