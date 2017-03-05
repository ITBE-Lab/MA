/* Keeps all configuration related info of the bio solution.
 */
#pragma once

#include <string>
#include <array>
#include <map>
#include <exception>
#include <iostream>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/expressions/attr_fwd.hpp>
#include <boost/log/expressions/attr.hpp>

#include <boost/phoenix.hpp>

#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>


class Configuration
{
public :
	boost::property_tree::ptree xPropertyTree;

	/* Delivers the folder for the storage of temporary data.
	 */
	boost::filesystem::path getTemporaryDataFolder( void )
	{
		return xPropertyTree.get<boost::filesystem::path>( "configure.global.prefix" )
			/= xPropertyTree.get<boost::filesystem::path>( "configure.global.temporary.folder" ); // name of the folder with taxonomy information
	} // method
	
	/* Delivers the full prefix for the NCBI database folder. (This folder keeps currently all sources as well)
	 */
	boost::filesystem::path getNCBIDatabasePrefix( void )
	{
		return
			   xPropertyTree.get<boost::filesystem::path>( "configure.global.prefix" )			// global prefix for data storage
			/= xPropertyTree.get<boost::filesystem::path>( "configure.NCBI.prefix" )			// NCBI section within the global data storage
			/= xPropertyTree.get<boost::filesystem::path>( "configure.NCBI.database.prefix" );	// database section within the NCBI section
	} // method

	boost::filesystem::path getNCBICacheFolder( void )
	{
		return
			   xPropertyTree.get<boost::filesystem::path>( "configure.global.prefix" )		// global prefix for data storage
			/= xPropertyTree.get<boost::filesystem::path>( "configure.NCBI.prefix" )		// NCBI section within the global data storage
			/= xPropertyTree.get<boost::filesystem::path>( "configure.NCBI.cache.folder" );	// cache folder within the NCBI folder
	} // method

	boost::filesystem::path getNCBIXMLFilesCacheFolder( void )
	{
		return getNCBICacheFolder()
			/= xPropertyTree.get<boost::filesystem::path>( "configure.NCBI.cache.xml.subfolder" ); // xml subfolder within the cache folder
	} // method

	boost::filesystem::path getNCBIFastaFilesCacheFolder( void )
	{
		return getNCBICacheFolder()
			/= xPropertyTree.get<boost::filesystem::path>( "configure.NCBI.cache.fasta.subfolder" ); // fasta subfolder within the cache folder
	} // method
	
	/* Delivers the fully qualified name of the database keeping all NCBI related data.
	 */
	boost::filesystem::path getNCBIDatabaseFileName( void )
	{
		return getNCBIDatabasePrefix()	// NCBI database prefix
			/= xPropertyTree.get<boost::filesystem::path>( "configure.NCBI.database.name" ); // name of the database
	} // method

	/* Delivers the fully qualified name of the folder that keeps the taxonomy sources
	 */
	boost::filesystem::path getTaxonomySourceFolder( void )
	{
		return getNCBIDatabasePrefix()	// NCBI database prefix
			/= xPropertyTree.get<boost::filesystem::path>( "configure.NCBI.database.taxonomy.prefix" ); // name of the folder with taxonomy information
	} // method

	boost::filesystem::path getAssembliesFolder( void )
	{
		return xPropertyTree.get<boost::filesystem::path>( "configure.global.prefix" )
			/= xPropertyTree.get<boost::filesystem::path>( "configure.NCBI.prefix" )
			/= xPropertyTree.get<boost::filesystem::path>( "configure.NCBI.assemblies.folder" );
	} // method

	/* Fully qualified name of the folder that holds all NCBI assembly sources.
	 * (Folder for keeping the .stat and .txt files for assemblies.)
	 */
	boost::filesystem::path getAssembliesSourcesFolder( void ) 
	{
		return getAssembliesFolder( )
			/= xPropertyTree.get<boost::filesystem::path>( "configure.NCBI.assemblies.sources.folder" );
	} // method

	boost::filesystem::path getAssembliesCacheFolder( void )
	{
		return xPropertyTree.get<boost::filesystem::path>( "configure.NCBI.assemblies.cache.folder" );
	} // method

	boost::filesystem::path getAssemblyPacksFolder( void )
	{
		return getAssembliesCacheFolder( )
			/= xPropertyTree.get<boost::filesystem::path>( "configure.NCBI.assemblies.cache.packs.subfolder" );
	} // method

	/* This folder is completely described in a single XML entry
	 */
	boost::filesystem::path getFM_IndicesFolder( void )
	{	return getAssembliesCacheFolder()
			/= xPropertyTree.get<boost::filesystem::path>( "configure.NCBI.assemblies.cache.fm-indices.subfolder" );
	} // method

	/* This folder contains preassembled FASTA collections for assemblies.
	 */
	boost::filesystem::path getAssemblyPreassembledFolder( void )
	{
		return getAssembliesFolder( )
			/= xPropertyTree.get<boost::filesystem::path>( "configure.NCBI.assemblies.preassembled.subfolder" );
	} // method

	boost::filesystem::path gene2accesion_path_local( void )
	{
		return getNCBIDatabasePrefix() /= "gene2accession";
	} // method

	boost::filesystem::path geneInfo_path_NCBI( void )
	{
		return getNCBIDatabasePrefix() /= "gene_info";
	} // method

	boost::filesystem::path geneNeighbors_path_local( void )
	{
		return getNCBIDatabasePrefix() /= "gene_neighbors";
	} // method

	boost::filesystem::path assemblySummaryGenbank_path_local( void )
	{
		return getNCBIDatabasePrefix() /= "assembly_summary_genbank.txt";
	} // method

	boost::filesystem::path taxonomyDumpNodesPath( void )
	{
		return getTaxonomySourceFolder() /= "nodes.dmp";
	} // method

	boost::filesystem::path taxonomyDumpNamesPath( void )
	{
		return getTaxonomySourceFolder() /= "names.dmp";
	} // method

	/* The creation of the map is quite inefficient, but we do this extremely seldom.
	 */
	std::map<std::string, std::string> getNCBIDatabaseSources()
	{
		std::map<std::string, std::string> xLocationMap;
		for ( auto &rxChild : xPropertyTree.get_child( "configure.NCBI.database.sources" ) )
		{
			xLocationMap[ rxChild.second.get<std::string>( "location" ) ] = rxChild.second.get<std::string>( "name" );
		} // for
		return xLocationMap;
	} // method

	/* This array controls the filter for trivial logging. true -> level is switeched on, false -> lebvel is switched off.
	 */
	std::array<bool, boost::log::trivial::fatal + 1> aLogLevelSwitches; //  = { { false, false, false, false, false, false } };

	/* Initializes the behavior of trivial logging.
	 * Defined in .cpp file.
	 */
	void vInitializeLogging( void );
	
	/* Constructor
	 */
	Configuration( const std::string &rsConfigurationAsXML ) :
		xPropertyTree(), // initialize the property tree
		aLogLevelSwitches({ { false, false, false, false, false, false } })
	{
		try 
		{	/* Feed a string stream with the configuration text.
			 */
			std::stringstream xStringStream;
			xStringStream << rsConfigurationAsXML;

			/* Parse the XML into the property tree.
			 */
			boost::property_tree::read_xml( xStringStream, xPropertyTree );
		} // try
		catch ( std::exception &rxException )
		{
			BOOST_LOG_TRIVIAL(fatal) << "Getting configuration failed due to reason: " << rxException.what();
			
			/* With an incorrect configuration is no value to continue.
			 */
			throw;
		} // catch

		/* We successfully parsed the configuration, so we are able to initialize logging.
		 */
		vInitializeLogging();

		/* Create several directories, specified in the configuration.
		 * FIX ME: 1) Maybe we should this do somewhere else. 2) If the directories does not exist, we get no exception at all.
		 */
		boost::filesystem::create_directories( getNCBIXMLFilesCacheFolder() );
		boost::filesystem::create_directories( getNCBIFastaFilesCacheFolder() );
		boost::filesystem::create_directories( getAssemblyPacksFolder() );
		boost::filesystem::create_directories( getFM_IndicesFolder() );
	} // constructor
}; // class

extern Configuration xGlobalConfiguration;

/* Delivers (creates, if necessary) a empty directory on the foundation of rsPath.
 * If the requested directory exists already, we create a fresh one by following a numbering scheme.
 * Maybe there is better location for for this function ...
 */
boost::filesystem::path sMakeDirectoryIfNotExists( const boost::filesystem::path &rsPath ); // prototype

/* Delivers the home directory of the current user for Linux and Windows.
 * Maybe there is better location for for this function ...
 */
boost::filesystem::path sGetHomeDirectory(); // prototype
