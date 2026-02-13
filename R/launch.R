# launch.R
## used to launch the shiny app in a web browser.


#' Launch shiny app
#'
#' @param verbose for debug use
#' @param ReductionKeyWords key words used for prepare Reduction options
#' @param SplitOptionMaxLevel max level cutoff for prepare Split options
#' @param MaxInputFileSize set the limited upload file size
#' @param dataset_dir primary server-side directory containing .rds/.qs2 datasets; app also scans ~/session_data/mounted-data-readonly
#' @param allow_browser_upload whether to show legacy browser upload input
#' @param workspace_dir base server directory used by "Save to Workspace" actions
#' @param allow_browser_download whether to keep browser download buttons (legacy mode)
#'
#' @import shiny
#' @return In-browser Shiny Application launch
#' @examples
#' if(interactive()){launchSeuratExplorer()}
#' @export
launchSeuratExplorer <- function(verbose = FALSE,
                                 ReductionKeyWords = c("umap","tsne"),
                                 SplitOptionMaxLevel = 12,
                                 MaxInputFileSize = 20*1024^3, # default 20GB
                                 dataset_dir = "data",
                                 allow_browser_upload = FALSE,
                                 workspace_dir = getwd(),
                                 allow_browser_download = FALSE
                                 ){
  options(SeuratExplorerVerbose = verbose)
  options(SeuratExplorerReductionKeyWords = ReductionKeyWords)
  options(SeuratExplorerSplitOptionMaxLevel = SplitOptionMaxLevel)
  options(shiny.maxRequestSize = MaxInputFileSize)
  options(SeuratExplorerDatasetDir = dataset_dir)
  options(SeuratExplorerAllowBrowserUpload = allow_browser_upload)
  options(SeuratExplorerWorkspaceDir = workspace_dir)
  options(SeuratExplorerAllowBrowserDownload = allow_browser_download)

  shinyApp(ui, server)
}
