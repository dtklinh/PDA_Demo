data <- psmelt(p.true.filt.s.unify)

custom_rarefy <- function(x,threshold){
  x <- x[x>0]
  sum(1-exp(lchoose(sum(x)-x,threshold) - lchoose(sum(x),threshold)))
}

data %>% 
  dplyr::group_by(Sample) %>% 
  dplyr::summarise(s_r = custom_rarefy(Abundance,1000), 
            s_ob = sum(Abundance >0), 
            n = sum(Abundance)) %>% 
  pivot_longer(cols = c(s_r,s_ob),values_to = "s",names_to = "type") %>% 
  ggplot(aes(x=n,y=s,color = type)) + scale_x_log10() +
  geom_point() +
  geom_smooth()

p.rare <- rarefy_even_depth(p.true.RmLowAbun.RmHighPrev.man,sample.size = 100,rngseed = 4)
beta_analysis_bray(p.rare)

bray_dist = phyloseq::distance(Phyloseq, method="bray")
p.adonis <- pairwise.adonis(x = bray_dist,factors = meta(Phyloseq)[,i])

p_hack <- function(Phyloseq,sample_size,seed,condition){
  p.rare <- rarefy_even_depth(Phyloseq,sample.size = sample_size,rngseed = seed)
  bray_dist = phyloseq::distance(p.rare, method="bray")
  p.adonis <- pairwise.adonis(x = bray_dist,factors = meta(p.rare)[,condition])
  return(p.adonis$p.adjusted)
}
vec <- seq(1,10000)
for (j in vec) {
  vec[j] <- p_hack(Phyloseq = p.true.RmLowAbun.RmHLPrev.man,sample_size = j,seed = 4,condition = "sample_side")
  print(j)
}
df_phack <- data.frame(x = seq(1,10000),
                       p = vec)
pp <- df_phack %>% ggplot(aes(y = p,x=x)) + geom_point()
pp
rm(j)
vec2 <- seq(10001,30000)
for (j in vec2) {
  vec2[j] <- p_hack(Phyloseq = p.true.RmLowAbun.RmHLPrev.man,sample_size = j,seed = 4,condition = "sample_side")
  print(j)
}
df_phack2 <- data.frame(x = seq(10001,30971),
                       p = vec2)
pp2 <- df_phack2 %>% ggplot(aes(y = p,x=x)) + geom_point()
pp2