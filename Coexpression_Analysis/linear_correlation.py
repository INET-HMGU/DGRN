import  matplotlib.pyplot   as  plt
import  numpy               as  np
import  pandas              as  pd
import  scipy.stats         as  stats

subset          =   np.zeros( ( 0 ) )
reference       =   np.zeros( ( 0 ) )

for species     in  [ "Ath", "Aly", "Esa" ]:
    
    TPM             =   pd.read_excel( "Log2FC.xlsx", sheet_name = species, index_col = 0, usecols = [ 1, 2, 3, 4 ] )
    interactions    =   pd.read_excel( "Interaction_list.xlsx", sheet_name = species + "_Interactions" )
    
    names           =   np.asarray( TPM.index,    dtype = str   )
    TPM             =   np.asarray( TPM,          dtype = float )
    interactions    =   np.asarray( interactions, dtype = str   )
    
    if species      ==   "Esa":
        names           =   np.asarray( names, dtype = "<U15" )
    
    indexA          =   np.where( [ names == entry for entry in np.unique( interactions[ :, 0 ] ) ] )[ 1 ]
    indexB          =   np.where( [ names == entry for entry in np.unique( interactions[ :, 1 ] ) ] )[ 1 ]
    
    A               =   TPM[ indexA ]
    B               =   TPM[ indexB ]
    
    A               =   np.asarray( [ A.T ** 0, A.T ** 1 ] ).T
    linalg          =   np.asarray( [ np.linalg.lstsq( a, B.T, rcond = None )[ 0 ].T for a in A ] )
    mse             =   np.mean( ( np.sum( A[ :, np.newaxis, :, : ] * linalg[ :, :, np.newaxis, : ], axis = 3 ) - B[ np.newaxis, :, : ] ) ** 2, axis = 2 )
    r2              =   1 - mse / np.var( B, axis = 1 )
    
    namesA          =   names[ indexA ]
    namesB          =   names[ indexB ]
    
    A               =   [ np.where( namesA == entry )[ 0 ][ 0 ] for entry in interactions[ :, 0 ] ]
    B               =   [ np.where( namesB == entry )[ 0 ][ 0 ] for entry in interactions[ :, 1 ] ]
    
    index           =   np.zeros( ( len( indexA ), len( indexB ) ), dtype = bool )
    index[ A, B ]   =   True
    
    subset          =   np.concatenate( ( subset,     r2[   index ] ) )
    reference       =   np.concatenate( ( reference,  r2[ ~ index ] ) )

subset          =   subset   [ np.isfinite( subset     ) ]
reference       =   reference[ np.isfinite( reference  ) ]

bins            =   np.arange( 0, 1.025, 0.025 )
hist_subset     =   np.histogram( subset,    bins = bins )[ 0 ]
hist_reference  =   np.histogram( reference, bins = bins )[ 0 ]

p_value         =   stats.mannwhitneyu( subset, reference ).pvalue



figure          =   plt.figure( figsize = ( 3.2, 3.2 ), layout = "constrained" )
ax              =   figure.add_subplot( 1, 1, 1 )

_               =   ax.boxplot( [ subset, reference ], positions = [ 1, 2 ], notch = True, widths = 0.2, whis = 0.05, showfliers = False )
_               =   ax.text( 1.5, 1.0, "p = %.2E" % p_value, horizontalalignment = "center" )

_               =   ax.set_xlim( [ 0.5, 2.5 ] )
_               =   ax.set_ylim( [ 0.0, 1.0 ] )
_               =   ax.set_xticks( ticks = [ 1, 2 ], labels = [ "Interacting\npairs", "Non-interacting\npairs" ] )
_               =   ax.set_yticks( ticks = [ 0, 0.2, 0.4, 0.6, 0.8, 1.0 ] )
_               =   ax.set_ylabel( "Coefficient of determination" )
_               =   ax.spines[ [ "top", "right" ] ].set_visible( False )

_               =   figure.savefig( "boxplot.pdf" )
_               =   figure.savefig( "boxplot.svg" )



figure          =   plt.figure( figsize = ( 3.2, 3.2 ), layout = "constrained" )
ax              =   figure.add_subplot( 1, 1, 1 )

_               =   ax.plot( ( bins[ : -1 ] + bins[ 1 : ] ) / 2, np.cumsum( hist_subset    ) / np.sum( hist_subset    ) )
_               =   ax.plot( ( bins[ : -1 ] + bins[ 1 : ] ) / 2, np.cumsum( hist_reference ) / np.sum( hist_reference ) )
_               =   ax.text( 0.5, 1.0, "p = %.2E" % p_value, horizontalalignment = "center" )

_               =   ax.set_xlim( [ 0.0, 1.0 ] )
_               =   ax.set_ylim( [ 0.0, 1.0 ] )
_               =   ax.set_xticks( ticks = [ 0, 0.2, 0.4, 0.6, 0.8, 1.0 ] )
_               =   ax.set_yticks( ticks = [ 0, 0.2, 0.4, 0.6, 0.8, 1.0 ] )
_               =   ax.set_xlabel( "Coefficient of determination" )
_               =   ax.set_ylabel( "Cumulative frequency" )
_               =   ax.legend( labels = [ "Interacting pairs", "Non-interacting pairs" ], ncol = 1, loc = 4, borderpad = 0, frameon = False )
_               =   ax.spines[ [ "top", "right" ] ].set_visible( False )

_               =   figure.savefig( "lineplot.pdf" )
_               =   figure.savefig( "lineplot.svg" )

#plt.show()

#p_value        =    stats.ttest_ind( subset[ np.isfinite( subset ) ], reference[ np.isfinite( reference ) ] ).pvalue

#print( "Ttest: " + str( p_value ) )

#p_value         =   stats.wilcoxon( subset[ np.isfinite( subset ) ], reference[ np.isfinite( reference ) ] ).pvalue

#print( "Wilcoxon: " + str( p_value ) )


#print( "mannwhitneyu: " + str( p_value ) )



#p_value        =    stats.ttest_ind( reference, subset, equal_var = True, alternative = "two-sided", nan_policy = "omit").pvalue
# 0.001276895665375408

# bins            =   np.arange( 0, 1.025, 0.025 )
# counts_1        =   np.histogram( subset,    bins = bins )[ 0 ]
# counts_2        =   np.histogram( reference, bins = bins )[ 0 ]

# figure          =   plt.figure( )
# ax_01           =   figure.add_subplot( 3, 2, 1 )
# ax_02           =   figure.add_subplot( 3, 2, 3 )
# ax_03           =   figure.add_subplot( 3, 2, 5 )

# ax_04           =   figure.add_subplot( 3, 2, 2 )
# ax_05           =   figure.add_subplot( 3, 2, 4 )
# ax_06           =   figure.add_subplot( 3, 2, 6 )

# _               =   ax_01.plot( ( bins[ : -1 ] + bins[ 1 : ] ) / 2, ( counts_1 ) / np.sum( counts_1 ) )
# _               =   ax_01.plot( ( bins[ : -1 ] + bins[ 1 : ] ) / 2, ( counts_2 ) / np.sum( counts_2 ) )

# _               =   ax_02.bar( ( bins[ : -1 ] + bins[ 1 : ] ) / 2 - 0.005, ( counts_1 ) / np.sum( counts_1 ), width = 0.0075 )
# _               =   ax_02.bar( ( bins[ : -1 ] + bins[ 1 : ] ) / 2 + 0.005, ( counts_2 ) / np.sum( counts_2 ), width = 0.0075 )

# _               =   ax_03.scatter( ( bins[ : -1 ] + bins[ 1 : ] ) / 2, ( counts_1 ) / np.sum( counts_1 ) )
# _               =   ax_03.scatter( ( bins[ : -1 ] + bins[ 1 : ] ) / 2, ( counts_2 ) / np.sum( counts_2 ) )

# _               =   ax_04.plot( ( bins[ : -1 ] + bins[ 1 : ] ) / 2, np.cumsum( counts_1 ) / np.sum( counts_1 ) )
# _               =   ax_04.plot( ( bins[ : -1 ] + bins[ 1 : ] ) / 2, np.cumsum( counts_2 ) / np.sum( counts_2 ) )

# _               =   ax_05.bar( ( bins[ : -1 ] + bins[ 1 : ] ) / 2 - 0.005, np.cumsum( counts_1 ) / np.sum( counts_1 ), width = 0.0075 )
# _               =   ax_05.bar( ( bins[ : -1 ] + bins[ 1 : ] ) / 2 + 0.005, np.cumsum( counts_2 ) / np.sum( counts_2 ), width = 0.0075 )

# _               =   ax_06.scatter( ( bins[ : -1 ] + bins[ 1 : ] ) / 2, np.cumsum( counts_1 ) / np.sum( counts_1 ) )
# _               =   ax_06.scatter( ( bins[ : -1 ] + bins[ 1 : ] ) / 2, np.cumsum( counts_2 ) / np.sum( counts_2 ) )

# _               =   ax_01.legend( labels = [ "+", "-" ] )
# _               =   ax_02.legend( labels = [ "+", "-" ] )
# _               =   ax_03.legend( labels = [ "+", "-" ] )

# _               =   plt.show( )