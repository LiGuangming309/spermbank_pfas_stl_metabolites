# devtools::install_github("hms-dbmi/UpSetR")

library(UpSetR)


movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), header=T, sep=";" )
mutations <- read.csv(system.file("extdata", "mutations.csv", package = "UpSetR"), header=T, sep = ",")

movies

upset(movies,attribute.plots=list(gridrows=60,plots=list(list(plot=scatter_plot, x="ReleaseDate", y="AvgRating"),
                                                         list(plot=scatter_plot, x="ReleaseDate", y="Watches"),list(plot=scatter_plot, x="Watches", y="AvgRating"),
                                                         list(plot=histogram, x="ReleaseDate")), ncols = 2))


mutations

upset(mutations, sets = c("PTEN", "TP53", "EGFR", "PIK3R1", "RB1"), sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on")


type(mutations)


names(mutations)


names(movies)

upset(movies, attribute.plots=list(gridrows = 100, ncols = 1, 
                                   plots = list(list(plot=histogram, x="AvgRating",queries=T),
                                                list(plot = scatter_plot, y = "AvgRating", x = "Watches", queries = T))), 
      sets = c("Action", "Adventure", "Children", "War", "Noir"),
      queries = list(list(query = intersects, params = list("War"), active = T),
                     list(query = intersects, params = list("Noir"))))




################################################################################
################################################################################

# install.packages("ggplot2")
# install.packages("ComplexUpset")
# install.packages("ggplot2movies")
# if(!require(devtools)) install.packages("devtools")
# devtools::install_github("krassowski/complex-upset")

library(ggplot2)
library(ComplexUpset)
# install.packages("ggplot2movies")
library(ggplot2movies)

ggplot2movies::movies

movies = as.data.frame(ggplot2movies::movies)
head(movies, 3)

movies$mpaa

genres = colnames(movies)[18:24]
genres


movies[genres] = movies[genres] == 1
t(head(movies[genres], 3))


movies[movies$mpaa == '', 'mpaa'] = NA
movies = na.omit(movies)




set_size = function(w, h, factor=1.5) {
  s = 1 * factor
  options(
    repr.plot.width=w * s,
    repr.plot.height=h * s,
    repr.plot.res=100 / factor,
    jupyter.plot_mimetypes='image/png',
    jupyter.plot_scale=1
  )
}



set_size(8, 3)
upset(movies, genres, name='genre', width_ratio=0.1)



set_size(8, 3)
(
  upset(movies, genres, name='genre', width_ratio=0.1, min_size=10, wrap=TRUE, set_sizes=FALSE)
  + ggtitle('Without empty groups (Short dropped)')
  +    # adding plots is possible thanks to patchwork
    upset(movies, genres, name='genre', width_ratio=0.1, min_size=10, keep_empty_groups=TRUE, wrap=TRUE, set_sizes=FALSE)
  + ggtitle('With empty groups')
)



set_size(8, 3)
upset(
  movies, genres, width_ratio=0.1,
  min_degree=3,
)


set_size(8, 3)
upset(
  movies, genres, width_ratio=0.1,
  n_intersections=15
)




abc_data = create_upset_abc_example()

abc_venn = (
  ggplot(arrange_venn(abc_data))
  + coord_fixed()
  + theme_void()
  + scale_color_venn_mix(abc_data)
)

(
  abc_venn
  + geom_venn_region(data=abc_data, alpha=0.05)
  + geom_point(aes(x=x, y=y, color=region), size=1)
  + geom_venn_circle(abc_data)
  + geom_venn_label_set(abc_data, aes(label=region))
  + geom_venn_label_region(
    abc_data, aes(label=size),
    outwards_adjust=1.75,
    position=position_nudge(y=0.2)
  )
  + scale_fill_venn_mix(abc_data, guide='none')
)


set_size(6, 6.5)
simple_venn = (
  abc_venn
  + geom_venn_region(data=abc_data, alpha=0.3)
  + geom_point(aes(x=x, y=y), size=0.75, alpha=0.3)
  + geom_venn_circle(abc_data)
  + geom_venn_label_set(abc_data, aes(label=region), outwards_adjust=2.55)
)
highlight = function(regions) scale_fill_venn_mix(
  abc_data, guide='none', highlight=regions, inactive_color='NA'
)

(
  (
    simple_venn + highlight(c('A-B')) + labs(title='Exclusive intersection of A and B')
    | simple_venn + highlight(c('A-B', 'A-B-C')) + labs(title='Inclusive intersection of A and B')
  ) /
    (
      simple_venn + highlight(c('A-B', 'A', 'B')) + labs(title='Exclusive union of A and B')
      | simple_venn + highlight(c('A-B', 'A-B-C', 'A', 'B', 'A-C', 'B-C')) + labs(title='Inclusive union of A and B')
    )
)





set_size(8, 4.5)
abc_upset = function(mode) upset(
  abc_data, c('A', 'B', 'C'), mode=mode, set_sizes=FALSE,
  encode_sets=FALSE,
  queries=list(upset_query(intersect=c('A', 'B'), color='orange')),
  base_annotations=list(
    'Size'=(
      intersection_size(
        mode=mode,
        mapping=aes(fill=exclusive_intersection),
        size=0,
        text=list(check_overlap=TRUE)
      ) + scale_fill_venn_mix(
        data=abc_data,
        guide='none',
        colors=c('A'='red', 'B'='blue', 'C'='green3')
      )
    )
  )
)

(
  (abc_upset('exclusive_intersection') | abc_upset('inclusive_intersection'))
  /
    (abc_upset('exclusive_union') | abc_upset('inclusive_union'))
)





set_size(8, 3)
upset(
  movies, genres,
  width_ratio=0.1,
  min_size=10,
  mode='inclusive_union',
  base_annotations=list('Size'=(intersection_size(counts=FALSE, mode='inclusive_union'))),
  intersections='all',
  max_degree=3
)


set_size(8, 8)

set.seed(0)   # keep the same jitter for identical plots

upset(
  movies,
  genres,
  annotations = list(
    # 1st method - passing list:
    'Length'=list(
      aes=aes(x=intersection, y=length),
      # provide a list if you wish to add several geoms
      geom=geom_boxplot(na.rm=TRUE)
    ),
    # 2nd method - using ggplot
    'Rating'=(
      # note that aes(x=intersection) is supplied by default and can be skipped
      ggplot(mapping=aes(y=rating))
      # checkout ggbeeswarm::geom_quasirandom for better results!
      + geom_jitter(aes(color=log10(votes)), na.rm=TRUE)
      + geom_violin(alpha=0.5, na.rm=TRUE)
    ),
    # 3rd method - using `upset_annotate` shorthand
    'Budget'=upset_annotate('budget', geom_boxplot(na.rm=TRUE))
  ),
  min_size=10,
  width_ratio=0.1
)



set_size(8, 5)

upset(
  movies,
  genres,
  annotations = list(
    'MPAA Rating'=(
      ggplot(mapping=aes(fill=mpaa))
      + geom_bar(stat='count', position='fill')
      + scale_y_continuous(labels=scales::percent_format())
      + scale_fill_manual(values=c(
        'R'='#E41A1C', 'PG'='#377EB8',
        'PG-13'='#4DAF4A', 'NC-17'='#FF7F00'
      ))
      + ylab('MPAA Rating')
    )
  ),
  width_ratio=0.1
)




set_size(8, 8)
set.seed(0)
upset(
  movies,
  genres,
  mode='inclusive_intersection',
  annotations = list(
    # if not specified, the mode will follow the mode set in `upset()` call (here: `inclusive_intersection`)
    'Length (inclusive intersection)'=(
      ggplot(mapping=aes(y=length))
      + geom_jitter(alpha=0.2, na.rm=TRUE)
    ),
    'Length (exclusive intersection)'=(
      ggplot(mapping=aes(y=length))
      + geom_jitter(alpha=0.2, na.rm=TRUE)
      + upset_mode('exclusive_intersection')
    ),
    'Length (inclusive union)'=(
      ggplot(mapping=aes(y=length))
      + geom_jitter(alpha=0.2, na.rm=TRUE)
      + upset_mode('inclusive_union')
    )
  ),
  min_size=10,
  width_ratio=0.1
)


upset_test(movies, genres)


chisq_from_formula = function(formula, data) {
  chisq.test(
    ftable(formula, data)
  )
}

anova_single = function(formula, data) {
  result = summary(aov(formula, data))
  list(
    p.value=result[[1]][['Pr(>F)']][[1]],
    method='Analysis of variance Pr(>F)',
    statistic=result[[1]][['F value']][[1]]
  )
}

custom_tests = list(
  mpaa=chisq_from_formula,
  budget=anova_single
)


head(upset_test(movies, genres, tests=custom_tests))



set_size(5, 3)

upset(
  movies, genres,
  min_size=10,
  width_ratio=0.3,
  set_sizes=(
    upset_set_size(
      geom=geom_bar(
        aes(fill=mpaa, x=group),
        width=0.8
      ),
      position='left'
    )
  ),
  # moves legends over the set sizes
  guides='over'
)


set_size(6, 4)
upset(
  movies,
  genres,
  min_size=10,
  width_ratio=0.2,
  stripes='white'
)









set_size(6, 4)
(
  upset(
    movies, genres, name='genre', width_ratio=0.1, min_size=10,
    stripes=c(alpha('grey90', 1), alpha('white', 1))
  )
  & theme(plot.background=element_rect(fill='white', color="white"))
)

set_size(8, 4)
upset(movies, genres, min_size=10, themes=list(default=theme()))

set_size(8, 8)

upset(
  movies,
  genres,
  annotations = list(
    'Length'=list(
      aes=aes(x=intersection, y=length),
      geom=geom_boxplot(na.rm=TRUE)
    ),
    'Rating'=list(
      aes=aes(x=intersection, y=rating),
      geom=list(
        geom_jitter(aes(color=log10(votes)), na.rm=TRUE),
        geom_violin(alpha=0.5, na.rm=TRUE)
      )
    )
  ),
  min_size=10,
  width_ratio=0.1,
  themes=modifyList(
    upset_themes,
    list(Rating=theme_void(), Length=theme())
  )
)


set_size(8, 4)

upset(
  movies, genres,
  base_annotations=list('Intersection size'=intersection_size(counts=FALSE)),
  min_size=100,
  width_ratio=0.1,
  themes=upset_modify_themes(
    list(
      'intersections_matrix'=theme(text=element_text(size=20)),
      'overall_sizes'=theme(axis.text.x=element_text(angle=90))
    )
  )
)



set_size(8, 6)

upset(
  movies, genres, name='genre', width_ratio=0.1, min_size=10,
  annotations = list(
    'Length'=list(
      aes=aes(x=intersection, y=length),
      geom=geom_boxplot(na.rm=TRUE)
    )
  ),
  queries=list(
    upset_query(
      intersect=c('Drama', 'Comedy'),
      color='red',
      fill='red',
      only_components=c('intersections_matrix', 'Intersection size')
    ),
    upset_query(
      set='Drama',
      fill='blue'
    ),
    upset_query(
      intersect=c('Romance', 'Comedy'),
      fill='yellow',
      only_components=c('Length')
    )
  )
)



























































