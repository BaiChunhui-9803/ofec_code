if (length(args) == 0) {
  stop(
    "Use: Rscript Data2GraphViz.r input_file output_folder \n  Example: Rscript Data2GraphViz.r ../Result/SpreadSpectrumRadarPollyPhase/data_filtered/data_SpreadSpectrumRadarPollyPhasen13_p1.256600.txt ../Result/SpreadSpectrumRadarPollyPhase/LONs/",
    call. = FALSE
  )
} else {
  # default output file
  file = args[1]
  output_folder = args[2]
}