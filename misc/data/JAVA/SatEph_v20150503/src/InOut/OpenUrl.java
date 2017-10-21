/*
 * OpenUrl.java
 *
 * Created on Utorok, 2008, april 15, 0:15
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 *
 * NOTES:  @author http://forum.java.sun.com/thread.jspa?threadID=5150750&messageID=9566398
 */

package InOut;

import javax.swing.JOptionPane;

/**
 *
 */
public class OpenUrl {
    
    /** Creates a new instance of OpenUrl */
    public OpenUrl() {
    }
    
    public void openUrl(String url){
     String osName = System.getProperty("os.name");
      try {
         if (osName.startsWith("Mac OS")) {
            //Class fileMgr = Class.forName("com.apple.eio.FileManager");
            //Method openURL = fileMgr.getDeclaredMethod("openURL",
             //  new Class[] {String.class});
            //openURL.invoke(null, new Object[] {url});
            }
         else if (osName.startsWith("Windows"))
            Runtime.getRuntime().exec("rundll32 url.dll,FileProtocolHandler " + url);
         else { //assume Unix or Linux
            String[] browsers = {
               "gnome-open", "firefox", "opera", "konqueror", "epiphany", "mozilla", "netscape" };
            String browser = null;
            for (int count = 0; count < browsers.length && browser == null; count++)
               if (Runtime.getRuntime().exec(
                     new String[] {"which", browsers[count]}).waitFor() == 0)
                  browser = browsers[count];
            if (browser == null)
               throw new Exception("Could not find web browser");
            else
               Runtime.getRuntime().exec(new String[] {browser, url});
            }
         }
      catch (Exception e) {
         JOptionPane.showMessageDialog(null, e.getLocalizedMessage());
         }
      }
    
}
