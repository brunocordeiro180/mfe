#' Clustering_cluster Meta-features
#'
#' Clustering_cluster measures extract information about validation index.
#'
#' @family meta-features
#' @param x A data.frame contained only the input attributes.
#' @param y A factor response vector with one label for each row/component of x.
#' @param features A list of features names or \code{"all"} to include all them.
#' @param summary A list of summarization functions or empty for all values. See
#'  \link{post.processing} method to more information. (Default: 
#'  \code{c("mean", "sd")})
#' @param transform A logical value indicating if the categorical attributes
#'  should be transformed. If \code{FALSE} they will be ignored. (Default: 
#'  \code{TRUE})
#' @param formula A formula to define the class column.
#' @param data A data.frame dataset contained the input attributes and class.
#'  The details section describes the valid values for this group.
#' @param ... Further arguments passed to the summarization functions.
#' @details
#'  The following features are allowed for this method:
#'  \describe{
#'    \item{"vdu"}{Calculate the Dunn Index.}
#'    \item{"vdb"}{Calculate the Davies and Bouldin Index.}
#'    \item{"int"}{Calculate the INT index.}
#'    \item{"sil"}{Calculate the mean silhouette value from data.}
#'    \item{"pb"}{Pearson Correlation between class matching and instance 
#'      distances.}
#'    \item{"ch"}{Calinski and Harabaz index.}
#'    \item{"nre"}{Normalized relative entropy.}
#'    \item{"sc"}{Mean of the number of examples per class.}
#'    \item{"xb"}{Ratio of overall deviation to cluster separation.}
#'  }
#' @return A list named by the requested meta-features.
#'
#' @references

#'  Bruno A. Pimentel, and Andre C. P. L. F. de Carvalho. A new data 
#'  characterization for selecting clustering_cluster algorithms using meta-learning. 
#'  Information Sciences, volume 477, pages 203 - 219, 2019.
#'
#' @examples
#' ## Extract all meta-features using formula
#' clustering_cluster(Species ~ ., iris)
#'
#' ## Extract some meta-features
#' clustering_cluster(iris[1:4], iris[5], c("vdu", "vdb", "sil"))
#'
#' ## Use another summarization function
#' clustering_cluster(Species ~ ., iris, summary=c("min", "median", "max"))
#' @export
clustering_cluster <- function(...) {
  UseMethod("clustering_cluster")
}

#' @rdname clustering_cluster
#' @export
clustering_cluster.default <- function(x, k=3, features="all",
                               summary=c("mean", "sd"),
                               transform=TRUE, ...) {
  if(!is.data.frame(x)) {
    stop("data argument must be a data.frame")
  }

  if(features[1] == "all") {
    features <- ls.clustering_cluster()
  }
  features <- match.arg(features, ls.clustering_cluster(), TRUE)
  colnames(x) <- make.names(colnames(x), unique=TRUE)

  if (length(summary) == 0) {
    summary <- "non.aggregated"
  }

  if(transform) {
    x <- binarize(x)
  } else {
    x <- x[sapply(x, is.numeric)]
  }

  x <- as.matrix(x)

  clusters <- kmeans(x, centers = k, nstart = 25)$cluster

  sapply(features, function(f) {
    fn <- paste("m", f, sep=".")
    measure <- eval(call(fn, x=x, y=clusters))
    post.processing(measure, summary, f %in% ls.clustering_cluster.multiples(), ...)
  }, simplify=FALSE)
}

#' @rdname clustering_cluster
#' @export
clustering_cluster.formula <- function(formula, data, features="all",
                                   summary=c("mean", "sd"),
                                   transform=TRUE, ...) {
  if(!inherits(formula, "formula")) {
    stop("method is only for formula datas")
  }

  if(!is.data.frame(data)) {
    stop("data argument must be a data.frame")
  }

  modFrame <- stats::model.frame(formula, data)
  attr(modFrame, "terms") <- NULL

  clustering_cluster.default(modFrame[-1], modFrame[1], features, summary, transform, 
    ...)
}

#' List the best clustering_cluster meta-features
#'
#' @return A list of best neighbor meta-features names.
#' @export
#'
#' @examples
#' ls.clustering_cluster()
ls.clustering_cluster <- function() {
  c("vdu", "vdb", "int", "sil", "pb", "ch", "nre", "sc", "xb", "knn_out", "bic", "aic", "c_index", "cm", "cn")
}

ls.clustering_cluster.multiples <- function() {
  c()
}

mc.vdu <- function(x, clusters) {
  aux <- clusterCrit::intCriteria(x, clusters, "Dunn")
  aux$dunn
}

mc.vdb <- function(x, clusters) {
  aux <- clusterCrit::intCriteria(x, clusters, "Davies_Bouldin")
  aux$davies_bouldin
}

mc.int <- function(x, clusters) {

  dfs <- ovo(x, factor(clusters))
  dst <- lapply(dfs, function(i) {
    dist(i$x)
  })

  aux <- mapply(function(dfs, dst) {
    inter(dfs, dst)
  }, dfs=dfs, dst=dst)

  c <- length(unique(y))
  aux <- sum(aux)/(c*(c-1)/2)
  return(aux)
}

mc.sil <- function(x, clusters) {
  aux <- clusterCrit::intCriteria(x, clusters, "Silhouette")
  aux$silhouette
}

mc.pb <- function(x, clusters) {
  aux <- clusterCrit::intCriteria(x, clusters, "Point_Biserial")
  aux$point_biserial
}

mc.ch <- function(x, clusters) {
  aux <- clusterCrit::intCriteria(x, clusters, "Calinski_Harabasz")
  aux$calinski_harabasz
}

mc.xb <- function(x, clusters) {
  aux <- clusterCrit::intCriteria(x, clusters, "Xie_Beni")
  aux$xie_beni
}

mc.nre <- function(x, clusters) {
  aux <- table(clusters)/length(clusters)
  -sum(aux * log2(aux))
}

mc.sc <- function(x, clusters) {
  mean(table(clusters))
}

mc.knn_out <- function(x, clusters){

  idxs <- sample(1:nrow(x),as.integer(0.7*nrow(x)), replace=FALSE)

  x_train <- x[idxs,]
  x_test <- x[-idxs, ]

  y_train <- clusters[idxs]
  y_test <- clusters[-idxs]

  ml <- class::knn(train= x_train, test = x_test, cl = y_train, k = 3)

  # confusion matrix
  confusion <- table(y_test,ml)

  # error
  length(y_test) - sum(y_test == ml)

}

mc.bic <- function(x, clusters){
  
  data <- as.data.frame(cbind(x, clusters))
  BIC(lm(data))
}

mc.aic <- function(x, clusters){
  
  data <- as.data.frame(cbind(x, clusters))
  AIC(lm(data))
}

mc.c_index <- function(x, clusters){
  
  aux <- clusterCrit::intCriteria(x, clusters, "C_index")
  aux$c_index
}

mc.cm <- function(x, clusters){

  data <- as.data.frame(cbind(x, clusters))

  dfs_by_label <- split(data, data$clusters)

  sum_distances <- c()

  for(df in dfs_by_label){
    centroid <- colMeans(df[, 1:ncol(x)])
    
    distances <- apply(df[,1:ncol(x)], 1, function(x){
      euc.dist(x, centroid)
    })
    
    sum_distances <- c(sum_distances, sum(distances))
  }

  mean(sum_distances)

}

mc.cn <- function(x, clusters){

  data <- as.data.frame(cbind(x, clusters))

  clValid::connectivity(data, data$clusters)

}
