SET(SRCs
  toolDisplay
  toolAlphaShape)

FOREACH(FILE ${SRCs})
  add_executable(${FILE} ${FILE})
  target_link_libraries( ${FILE} ${DGTAL_LIBRARIES} )
ENDFOREACH(FILE)
