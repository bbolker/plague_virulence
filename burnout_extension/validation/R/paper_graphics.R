# Common final-figure export: native LaTeX text, with plotmath converted before
# drawing so mathematical labels use the manuscript's Computer Modern fonts.
# Set TIKZ_XELATEX / TIKZ_TECTONIC when the compilers are not on PATH.
paper_math <- function(x) {
  if (is.expression(x)) return(vapply(as.list(x), function(z) paste0('$', paper_math(z), '$'), ''))
  if (is.numeric(x)) return(as.character(x))
  if (is.character(x)) return(x)
  if (is.symbol(x)) {
    s <- as.character(x)
    if (s %in% c('rho','theta','delta','gamma')) return(paste0('\\',s))
    return(s)
  }
  op <- as.character(x[[1]])
  args <- as.list(x)[-1]
  a <- vapply(args,paper_math,'')
  switch(op, '/' = paste0('\\frac{',a[1],'}{',a[2],'}'), '-'  = if(length(a)==1) paste0('-',a) else paste(a,collapse='-'), '[' = paste0(a[1],'_{',if(a[2] %in% c('in','BL')) paste0('\\mathrm{',a[2],'}') else a[2],'}'),
    '^' = paste0(a[1],'^{',a[2],'}'),
    '==' = paste(a,collapse=' = '), '<' = paste(a,collapse=' < '),
    'italic' = paste0('\\mathit{',a[1],'}'),
    'sqrt' = paste0('\\sqrt{',a[1],'}'),
    'paste' = paste0(a,collapse=''), '~' = paste(a,collapse='\\;'),
    '*' = paste0(a,collapse=''), paste(a,collapse=' '))
}
paper_text <- function(x) {
  if (is.call(x)) x <- as.expression(list(x))
  if (is.list(x)) return(unlist(lapply(x,paper_text),use.names=FALSE))
  if (is.expression(x)) x <- paper_math(x)
  x <- gsub('\u211b\u2080', '$\\mathcal{R}_0$', x, fixed=TRUE)
  x <- gsub('\u211b', '\\mathcal{R}', x, fixed=TRUE)
  x <- gsub('\u2264', '$\\leq$', x, fixed=TRUE)
  x <- gsub('\u2265', '$\\geq$', x, fixed=TRUE)
  x <- gsub('\u2212', '-', x, fixed=TRUE)
  x
}
paper_grob <- function(g) {
  if (inherits(g,'text')) {
    g$label <- paper_text(g$label)
    g$gp$fontfamily <- ''
  }
  if (!is.null(g$grobs)) g$grobs <- lapply(g$grobs,paper_grob)
  if (!is.null(g$children)) g$children <- do.call(grid::gList,lapply(g$children,paper_grob))
  g
}
paper_tikz_open <- function(file,width,height,pointsize=11) {
  compiler <- Sys.getenv('TIKZ_XELATEX',Sys.which('xelatex'))
  if (!nzchar(compiler)) stop('Set TIKZ_XELATEX to a XeLaTeX compiler.')
  options(tikzLatex=compiler,tikzXelatex=compiler,tikzDefaultEngine='xetex',
    tikzDocumentDeclaration=sprintf('\\documentclass[%dpt]{article}',pointsize),
    tikzMetricsDictionary=file.path(normalizePath('validation/figures/paper'),'.tikz-metrics'))
  tikzDevice::tikz(file,width=width,height=height,standAlone=TRUE,
    sanitize=FALSE,engine='xetex',pointsize=pointsize,
    packages=c('\\usepackage{tikz}','\\usepackage{amsmath,amssymb}',
               '\\usetikzlibrary{calc}',
               '\\usepackage[active,tightpage]{preview}',
               '\\PreviewEnvironment{tikzpicture}',
               '\\setlength{\\PreviewBorder}{0pt}'))
}
paper_compile <- function(tex) {
  tectonic <- Sys.getenv('TIKZ_TECTONIC',Sys.which('tectonic'))
  if (nzchar(tectonic)) {
    status <- system2(tectonic,c('--keep-logs','--outdir',shQuote(dirname(tex)),shQuote(tex)))
  } else {
    compiler <- Sys.getenv('TIKZ_XELATEX',Sys.which('xelatex'))
    status <- system2(compiler,c('-interaction=nonstopmode','-halt-on-error',
      '-output-directory',shQuote(dirname(tex)),shQuote(tex)))
  }
  if (status != 0) stop('Figure compilation failed: ',tex)
}
paper_save_plot <- function(p,name,width=6.8,height=5.5,landscape=FALSE,appendix=FALSE) {
  size <- if(landscape)17 else if(appendix)13 else 14
  p <- p + ggplot2::theme(text=ggplot2::element_text(family='',size=size),
    axis.text=ggplot2::element_text(size=if(appendix)10 else size-2),
    axis.title=ggplot2::element_text(size=size),
    strip.text=ggplot2::element_text(size=if(appendix)11 else size-1),
    legend.text=ggplot2::element_text(size=size-2),
    legend.title=ggplot2::element_text(size=size-1),
    axis.title.x=ggplot2::element_text(family=''),
    axis.title.y=ggplot2::element_text(family='',angle=90)) +
    ggplot2::labs(x='$\\mathcal{R}_0$')
  if (landscape) {
    # Regular grid in the transformed coordinates: rasterize only the colour field.
    p$layers[[1]] <- ggplot2::geom_raster(ggplot2::aes(R0,rho),interpolate=FALSE)
    p <- p + ggplot2::labs(y='Replenishment rate ($\\rho$)') +
      ggplot2::guides(fill=ggplot2::guide_colourbar(display='raster',order=1)) +
      ggplot2::theme(panel.spacing=grid::unit(.8,'lines'),strip.text.y.right=ggplot2::element_text(angle=90)) +
      ggplot2::scale_y_log10(breaks=c(.001,.005,.02,.05), labels=c('0.001','0.005','0.02','0.05'),expand=c(0,0)) +
      scale_x_log10_shifted(breaks=c(1.1,2,10),expand=c(0,0))
  }
  for (i in seq_along(p$scales$scales)) {
    s <- p$scales$scales[[i]]
    if (is.expression(s$name) && identical(as.character(s$name),'theta'))
      p$scales$scales[[i]]$name <- 'Recruitment shape ($\\theta$)'
  }
  # Translate plotmath labels before layout; native LaTeX math avoids
  # passing R's internal symbol-font codes to the TikZ text metrics.
  old_labeller <- p$facet$params$labeller
  if (!is.null(old_labeller)) p$facet$params$labeller <- function(labels) {
    lapply(old_labeller(labels), paper_text)
  }
  for (i in seq_along(p$scales$scales)) {
    s <- p$scales$scales[[i]]
    if (is.expression(s$name)) p$scales$scales[[i]]$name <- paper_text(s$name)
    if (is.expression(s$labels) || is.character(s$labels)) p$scales$scales[[i]]$labels <- paper_text(s$labels)
    if (is.function(s$labels)) {
      original_labels <- s$labels
      p$scales$scales[[i]]$labels <- local({ f <- original_labels; function(x) paper_text(f(x)) })
    }
  }
  if (grepl('unconditional_persistence_validation',name))
    p <- p + ggplot2::guides(colour=ggplot2::guide_legend(nrow=2))
  tex <- file.path('validation/figures/paper',paste0(name,'.tex'))
  paper_tikz_open(tex,width,height)
  grid::grid.draw(paper_grob(ggplot2::ggplotGrob(p)))
  grDevices::dev.off()
  paper_compile(tex)
}
