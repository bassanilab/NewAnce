package newance.mzjava.ms.spectrum;

import java.net.URI;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;

/**
 * @author fnikitin
 *         Date: 10/8/12
 */
public class URIBuilder {

    public static final URI UNDEFINED_URI = URI.create("/dev/null");

    private final String protocol = "software";
    private final String nameSpace;
    private final String appName;
    private String version;

    private final List<String> args;
    private final Properties params;

    public URIBuilder(String nameSpace, String appName) {

        this.nameSpace = nameSpace;
        this.appName = appName;

        this.args = new ArrayList<String>(10);
        this.params = new Properties();
    }

    public URIBuilder version(String version) {

        this.version = version;
        return this;
    }

    public URIBuilder arg(int index, String value) {

        args.add(index, value);
        return this;
    }

    public URIBuilder param(String key, String value) {

        params.put(key, value);
        return this;
    }

    public URI build() {

        if (nameSpace == null || nameSpace.length() == 0) {

            throw new IllegalStateException("name space is missing!");
        }

        if (appName == null || appName.length() == 0) {

            throw new IllegalStateException("software name is missing!");
        }

        StringBuilder sb = new StringBuilder();

        sb.append(protocol);
        sb.append("://");
        sb.append(nameSpace);
        sb.append("/");
        sb.append(appName);
        if (version != null) {

            sb.append(":v").append(version);
        }

        if (!args.isEmpty() || !params.isEmpty()) {

            sb.append("?");

            for (int i = 0; i < args.size(); i++) {

                sb.append("arg");
                sb.append(i);
                sb.append("=");
                sb.append(args.get(i));
                sb.append("&");
            }

            for (String key : params.stringPropertyNames()) {

                sb.append(key);
                sb.append("=");
                sb.append(params.getProperty(key));
                sb.append("&");
            }
            sb.delete(sb.length() - 1, sb.length());
        }

        try {

            return new URI(sb.toString());
        } catch (URISyntaxException e) {

            throw new IllegalStateException(e);
        }
    }
}
