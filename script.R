library(network)		# network data storage
library(sna)			# network analysis routines
library(latticeExtra)	# for nicer convergence & gof plots
library(ergm)			# fitting & evaluating ERGMs
library(igraph)


library(qgraph)
library(ggplot2)
library(glue)
# pre-processing
load("Pair3N.RData")
friendship <- network(fri)
G_fri <- graph_from_adjacency_matrix(fri) # directed 
fri_stat <- centrality(G_fri)
# add attributes to the nodes
friendship %v% 'gender' <- fem
friendship %v% 'gender.color' <- c('lightblue', 'pink')[fem + 1]
sch_normal <- (sch - mean(sch)) / sd(sch) + 4.5
friendship %v% 'school' <- sch

# set up parameters for the graph
set.seed(123)
pos <- plot(friendship)


# visualization of the network
plot(friendship, vertex.cex= sch_normal, vertex.col= 'gender.color',
     coord = pos)

# histogram - distribution of degree

df_out <- data.frame(fri_stat$OutDegree) 
ggplot(df_out, aes(x = fri_stat.OutDegree)) +
  geom_histogram(binwidth = 1, fill = "lightblue") +
  theme_light() +
  xlab('Out-Degree') +
  ylab('Count') +
  xlim(c(0, 10)) +
  ylim(c(0, 8))

df_in <- data.frame(fri_stat$InDegree) 
ggplot(df_in, aes(x = fri_stat.InDegree)) +
  geom_histogram(binwidth = 1, fill = "lightblue") +
  theme_light() +
  xlab('In-Degree') +
  ylab('Count') +
  xlim(c(0, 10)) +
  ylim(c(0, 8))

# transitivity 
# print(glue('Transitivity in this network is {round(gden(friendship), 2)}'))


library(texreg) # for showing ascII tables
library(kableExtra)
library(gtsummary)
library(gt)
# Q5
# full model
model_full <- ergm(friendship ~ edges + mutual + 
                     nodematch('gender') + gwesp(0.7,fixed=TRUE) +
                     twopath + gwidegree(0.7, fixed=TRUE))
# no closure tendencies 
model_noc <- ergm(friendship ~ edges + mutual + 
                    nodematch('gender') + 
                    gwidegree(0.7, fixed=TRUE))

# no Matthew effect 
model_matthew <- ergm(friendship ~ edges + mutual + nodematch('gender')+ gwesp(0.7,fixed=TRUE) + twopath)

# null model
model_nul <- ergm(friendship ~ edges + mutual + 
                    nodematch('gender'))

models <- list(model_full, model_noc, model_matthew, model_nul)
# makeing a table
# screenreg(models)
new_names <- c('edges' = 'edges',
               'mutual' = 'mutual',
               'nodematch.gender' = 'homophily*',
               'gwesp.OTP.fixed.0.7' = 'closure tend.*',
               'twopath' = 'two path',
               'gwideg.fixed.0.7' = 'preferential tend*')



library(texreg) # for showing ascII tables
library(modelsummary)
library(gt)

# change the columns names
new_names <- c('edges' = 'edges',
               'mutual' = 'mutual',
               'nodematch.gender' = 'homophily*',
               'gwesp.OTP.fixed.0.7' = 'closure tend.*',
               'twopath' = 'two path',
               'gwideg.fixed.0.7' = 'preferential tend*')
# make the table
modelsummary(models, stars = c('*' = .05, '**' = 0.01, 
                               '***' = 0.001),
             fmt = 2, exponentiate = TRUE,
             statistic = 'conf.int', 
             conf_level = .95, coef_map = new_names,
             output = 'gt') %>% 
  tab_spanner(label = "Full",columns = '(1)') %>%
  tab_spanner(label = md('No Clo<sup>1</sup>'),
              columns = '(2)') %>% 
  tab_spanner(label = md("No Matt<sup>2</sup>"), 
              columns = '(3)') %>% 
  tab_spanner(label = "Null",columns = '(4)') %>% 
  tab_header(title = 'Different Model Specification') %>% 
  tab_source_note(source_note = md('<sup>1</sup> No transitive closure tendencies ')) %>% 
  tab_source_note(source_note = md('<sup>2</sup> No Matthew effect') ) %>% 
  tab_source_note(source_note = md('<sup>*</sup> Refers to *Gender* homophily.')) %>% 
  tab_source_note(source_note = md('<sup>*</sup> Closure tendency ')) %>% 
  tab_source_note(source_note = md('<sup>*</sup> Preferential attachment tendencies')) %>% 
  tab_source_note(source_note = md("*Figures based on authors' calculation*")) %>% 
  tab_source_note(
    source_note = md('**Reference**: *Knecht, A. B. (2008). Friendship selection and friends influence. Dynamics of networks and actor attributes in early adolescence. PhD dissertation. Utrecht University*'
    ))


# empirical 
# par(mfrow = c(1, 1), bg = 'cornsilk')
plot(friendship, vertex.cex= sch_normal, vertex.col= 'gender.color',
     coord = pos)

# simulations
# full model
# par(mfrow = c(1, 1), bg = 'white')
sims_full <- simulate(model_full, nsim = 100, seed = 456)
plot(sims_full[[1]], vertex.cex= sch_normal, 
     vertex.col= 'gender.color', coord = pos)

# no closure
sims_noc <- simulate(model_noc, nsim = 100, seed = 456)
plot(sims_noc[[1]],  vertex.cex= sch_normal, 
     vertex.col= 'gender.color', coord = pos)

# no Matthew effect
sims_matthew <- simulate(model_matthew, nsim = 100, seed =456)
plot(sims_matthew[[1]],  vertex.cex= sch_normal, 
     vertex.col= 'gender.color', coord = pos)

# null model
sims_nul <- simulate(model_nul, nsim = 100, seed = 456)
plot(sims_nul[[1]],  vertex.cex= sch_normal, 
     vertex.col= 'gender.color', coord = pos)


# density - empirical value vs. different models 
par(mar = c(4, 4, 3, 9))
plot(density(gden(sims_full)), col = 'deepskyblue', lwd = 2, 
     main = '') 
lines(density(gden(sims_noc)), col = 'pink', lwd = 2)
lines(density(gden(sims_matthew)), col = 'darkorchid1', lwd = 2)
lines(density(gden(sims_nul)), col = 'aquamarine2', lwd = 2)
abline(v = gden(friendship), lty = "dashed", lwd = 1.5)
# legend('topright', inset = c(-0.57, 0) ,legend = c('Full', 'No Closure', 'No Matthew',
#                               'Null'), col = c('deepskyblue', 'pink','darkorchid1', 'aquamarine2'), lwd = 2, xpd = TRUE)
# transitivity 
par(mar = c(4, 4, 3, 9))
plot(density(gtrans(sims_full)), ylim = c(0, 8), 
     col = 'deepskyblue', lwd = 2, main = '')
lines(density(gtrans(sims_noc)), col = 'pink', lwd = 2)
lines(density(gtrans(sims_matthew)), col = 'darkorchid1', lwd = 2)
lines(density(gtrans(sims_nul)),col = 'aquamarine2', lwd = 2)
abline(v = gtrans(friendship), lty = "dashed", lwd = 1.5)
legend('topright', inset = c(-0.57, 0) ,legend = c('Full', 'No Closure', 'No Matthew',
                                                   'Null'), col = c('deepskyblue', 'pink','darkorchid1', 'aquamarine2'), lwd = 2, xpd = TRUE)



# goodness of model fit
# this is 
fit_no_matthew <- gof(model_matthew, 
                      control = control.gof.ergm(seed = 789, nsim=200))
plot(fit_no_matthew)

