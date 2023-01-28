def get_idf_file():
    # ------------------------- a) import library ------------------------ #
    from eppy.modeleditor import IDF
    # ----------------------------- b) names ----------------------------- #
    epfdir = '/home/fernando/software/EnergyPlus-22-2-0/'
    idf_file="1ZoneUncontrolled_win_1.idf"
    # -------------------------- c) set idd file ------------------------- #
    iddfile = epfdir + 'Energy+.idd'
    IDF.setiddname(iddfile)
    # ------------------------- d) get zone name ------------------------- #
    idf = IDF(epfdir+'ExampleFiles/'+idf_file)
    return idf